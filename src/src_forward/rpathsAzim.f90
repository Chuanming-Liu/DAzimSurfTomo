!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates ray path geometries for each
! source-receiver combination. It will also compute
! Frechet derivatives using these ray paths if required.
!input:
!scx, scz:  current source coordinates, in rad
!surfrcx, surfrcz= station coordinates, in rad
!writepath= path
!output:
!fmd=Frechet derivative matrix, 2D, i th data (fixed period, source, receiver)
! modified by Chuanming Liu, 2016, Oct.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE rpaths(wrgf,csid,cfd,scx,scz)
SUBROUTINE rpathsAzim(scx,scz,fdm,fdmc,fdms,surfrcx,surfrcz,writepath,Tperiod)
USE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER :: nopath=0
INTEGER :: i,j,k,l,m,n,ipx,ipz,ipxr,ipzr,nrp,sw
!fang!INTEGER :: wrgf,cfd,csid,ipxo,ipzo,isx,isz
INTEGER :: ipxo,ipzo,isx,isz
INTEGER :: ivx,ivz,ivxo,ivzo,nhp,maxrp
INTEGER :: ivxt,ivzt,ipxt,ipzt,isum,igref
INTEGER, DIMENSION (4) :: chp
REAL(KIND=i5) :: rayx,rayz
REAL(KIND=i10) :: dpl,rd1,rd2,xi,zi,vel,velo
REAL(KIND=i10) :: v,w,rigz,rigx,dinc,scx,scz
REAL(KIND=i10) :: dtx,dtz,drx,drz,produ,sred
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rgx,rgz
!fang!REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: fdm
REAL(KIND=i10), DIMENSION (4) :: vrat,vi,wi,vio,wio
!fang!------------------------------------------------
REAL fdm(0:nvz+1,0:nvx+1)
REAL(KIND=i10) surfrcx,surfrcz
LOGICAL:: writepath
!Liu!-------------------------------------------------
REAL fdmc(0:nvz+1,0:nvx+1)
REAL fdms(0:nvz+1,0:nvx+1)
REAL(KIND=i10) :: rdc1,rdc2
REAL(KIND=i10) :: rgx1,rgz1,rgx2,rgz2,delta,az,baz
REAL(KIND=i10) :: rgpsi
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rgaz
REAL*8 Tperiod
REAL(KIND=i10) :: x,z,goxb,gozb
INTEGER jj,kk
real,parameter::ftol=1e-4

!fang!------------------------------------------------
!
! ipx,ipz = Coordinates of cell containing current point, propagation coordinate (in grid)
! ipxr,ipzr = Same as ipx,apz except for refined grid
! ipxo,ipzo = Coordinates of previous point
! rgx,rgz = (x,z) coordinates of ray geometry, from station to source, in rad
! rgpsi = psi angle, azimuth from rgx(j) to rgx(j+1), in rad, segment azimuthal angle--rgaz in degree
! ivx,ivz = Coordinates of B-spline vertex containing current point, inversion grid
! ivxo,ivzo = Coordinates of previous point
! maxrp = maximum number of ray points, nnx, nnz grid number of propagation grid
! nrp = number of points to describe ray
! dpl = incremental path length of ray
! xi,zi = edge of model coordinates
! dtx,dtz = components of gradT
! wrgf = Write out raypaths? (<0=all,0=no,>0=souce id)
! cfd = calculate Frechet derivatives? (0=no,1=yes)
! csid = current source id
! fdm = Frechet derivative matrix,
! fdmc = Frechet derivative matrix, cos2pi part
! fdms = Frechet derivative matrix, sin2pi part
! nhp = Number of ray segment-B-spline cell hit points, Number of ray segments
! vrat = length ratio of ray sub-segment
! chp = pointer to incremental change in x or z cell
! drx,drz = distance from reference node of cell
! produ = variable for trilinear interpolation
! vel = velocity at current point, on ray path hit points
! velo = velocity at previous point
! v,w = local variables of x,z
! vi,wi = B-spline basis functions at current point
! vio,wio = vi,wi for previous point
! ivxt,ivzt = temporary ivr,ivx,ivz values
! rigx,rigz = end point of sub-segment of ray path
! ipxt,ipzt = temporary ipx,ipz values
! dinc = path length of ray sub-segment
! rayr,rayx,rayz = ray path coordinates in single precision
! isx,isz = current source cell location
! scx,scz = current source coordinates, in rad
! surfrcx,surfrcz= current receiver coordinates, in rad
! sred = source to ray endpoint distance, signal to judge exit the cycle.
! igref = ray endpoint lies in refined grid? (0=no,1=yes)
! nopath = switch to indicate that no path is present
! asgr = Apply source grid refinement (0=no,1=yes), set asgr =1
!
!-------------------------globalp summary--------------------------------------------------------------------
! Travel time field
! ttn(i,j) = traveltime field on the refined grid of nodes
! ttnr(i,j) = ttn for refined grid

! inversion grid
! dvx,dvz = B-spline vertex separation; inversion grid interval (in rad)
! dvxd,dvzd = dvx,dvz in degrees, inversion grid interval (degree)
! gox,goz = Origin of grid (theta,phi) (rad), origin of inversion grid (exchange in cal goxr, but retrieved soon)

! propagation grid
! nnx,nnz = Number of nodes of grid in x and z, number of propagation grid.
! gdx,gdz = grid dicing in x and z, number , set as gdx=gdz=5   grid point interval of propagation grid (in point)
! dnx,dnz = Node separation interval of grid in  x and z, in rad; dnx=dvx(velocity node)/ gdx,    dnx: grid interval of propagation grid (in rad)
! dnzd,dnzd = dnx,dnz in degrees,  propagation grid interval (degree)

! refined source grid
!  nnxr,nnzr = Number of nodes of refined grid in x and z
! goxr, gozr = Origin of refined grid (theta,phi) (in rad), refined source rectangular start point
! dnxr,dnzr = Node separation of refined grid in x and z (rad), dnxr=dvx/(gdx*sgdl)
! veln(i,j) = velocity values on a refined grid of nodes
! vnl,vnr,vnb,vnt = Bounds of refined grid
! nrnx,nrnz = Number of nodes in x and z for refined grid
! gorx,gorz = Grid origin of refined grid
!--------------------------------------------------------------------------------------------------------------------------------------------
! Allocate memory to arrays for storing ray path geometry
!
maxrp=nnx*nnz
ALLOCATE(rgx(maxrp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgx'
ENDIF
ALLOCATE(rgz(maxrp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgz'
ENDIF
ALLOCATE(rgaz(maxrp+1), STAT=checkstat)
rgx=0
rgz=0
rgaz=0
!
! Allocate memory to partial derivative array
!
!fang!IF(cfd.EQ.1)THEN
!fang!   ALLOCATE(fdm(0:nvz+1,0:nvx+1), STAT=checkstat)
!fang!   IF(checkstat > 0)THEN
!fang!      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL fdm'
!fang!   ENDIF
!fang!ENDIF
!
! Locate current source cell
!
IF(asgr.EQ.1)THEN
   isx=INT((scx-goxr)/dnxr)+1
   isz=INT((scz-gozr)/dnzr)+1
ELSE
   isx=INT((scx-gox)/dnx)+1
   isz=INT((scz-goz)/dnz)+1
ENDIF
!
! Set ray incremental path length equal to half width
! of cell
!
  dpl=dnx*earth
  rd1=dnz*earth*SIN(gox)
  IF(rd1.LT.dpl)dpl=rd1
  rd1=dnz*earth*SIN(gox+(nnx-1)*dnx)
  IF(rd1.LT.dpl)dpl=rd1
  dpl=0.5*dpl
!
! Loop through all the receivers
!
!fang!DO i=1,nrc
!
!  If path does not exist, then cycle the loop
!
fdm=0
! Liu
fdmc=0
fdms=0
!fang!   IF(cfd.EQ.1)THEN
!fang!      fdm=0.0
!fang!   ENDIF
!fang!   IF(srs(i,csid).EQ.0)THEN
!fang!      IF(wrgf.EQ.csid.OR.wrgf.LT.0)THEN
!fang!         WRITE(40)nopath
!fang!      ENDIF
!fang!      IF(cfd.EQ.1)THEN
!fang!         WRITE(50)nopath
!fang!      ENDIF
!fang!      CYCLE
!fang!   ENDIF
!------------------------------- Step1.---------------------------!
!  The first step is to locate the receiver in the grid.
!
   ipx=INT((surfrcx-gox)/dnx)+1
   ipz=INT((surfrcz-goz)/dnz)+1
   sw=0
   IF(ipx.lt.1.or.ipx.ge.nnx)sw=1
   IF(ipz.lt.1.or.ipz.ge.nnz)sw=1
   IF(sw.eq.1)then
      surfrcx=90.0-surfrcx*180.0/pi
      surfrcz=surfrcz*180.0/pi

         goxb=gox
         gozb=goz
         gox=90.0-goxb*180.0/pi
         goz=gozb*180.0/pi

      WRITE(6,*)"  Receiver lies outside model (lat,long)= ",surfrcx,surfrcz
     WRITE(6,*)"Boundary of model is (lat,lon) NE boundary= ",gox,goz
         gox=90.0-(goxb+(nnx-1)*dnx)*180.0/pi
         goz=(gozb+(nnz-1)*dnx)*180.0/pi
         WRITE(6,*)"Boundary of model is (lat,lon) SW boundary= ",gox,goz
      WRITE(6,*)"Index (lat,long)=",ipx,ipz
      WRITE(6,*)"Index limit (lat,long)=",nnx,nnz




         WRITE(6,*)"TERMINATING PROGRAM!!! in subroutine: rpathsAzim"

      STOP
   ENDIF
   IF(ipx.eq.nnx)ipx=ipx-1
   IF(ipz.eq.nnz)ipz=ipz-1
!
!  First point of the ray path is the receiver
!
   rgx(1)=surfrcx
   rgz(1)=surfrcz
!
!  Test to see if receiver is in source neighbourhood
!
   sred=((scx-rgx(1))*earth)**2
   sred=sred+((scz-rgz(1))*earth*SIN(rgx(1)))**2
   sred=SQRT(sred)
   IF(sred.LT.2.0*dpl)THEN
      rgx(2)=scx
      rgz(2)=scz
      nrp=2
      sw=1
   ENDIF
!
!  If required, see if receiver lies within refined grid
!
   IF(asgr.EQ.1)THEN
      ipxr=INT((surfrcx-goxr)/dnxr)+1
      ipzr=INT((surfrcz-gozr)/dnzr)+1
      igref=1
      IF(ipxr.LT.1.OR.ipxr.GE.nnxr)igref=0
      IF(ipzr.LT.1.OR.ipzr.GE.nnzr)igref=0
      IF(igref.EQ.1)THEN
         IF(nstsr(ipzr,ipxr).NE.0.OR.nstsr(ipzr+1,ipxr).NE.0)igref=0
         IF(nstsr(ipzr,ipxr+1).NE.0.OR.nstsr(ipzr+1,ipxr+1).NE.0)igref=0
      ENDIF
   ELSE
      igref=0
   ENDIF
!---------------------Step 2.-----------------------------------!
!  Due to the method for calculating traveltime gradient, if the
!  the ray end point lies in the source cell, then we are also done.
!
   IF(sw.EQ.0)THEN
      IF(asgr.EQ.1)THEN
         IF(igref.EQ.1)THEN
            IF(ipxr.EQ.isx)THEN
               IF(ipzr.EQ.isz)THEN
                  rgx(2)=scx
                  rgz(2)=scz
                  nrp=2
                  sw=1
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF(ipx.EQ.isx)THEN
            IF(ipz.EQ.isz)THEN
               rgx(2)=scx
               rgz(2)=scz
               nrp=2
               sw=1
            ENDIF
         ENDIF
      ENDIF
   ENDIF
!----------------------Step3.-------------------------------!
!  Now trace ray from receiver to "source"
!



   DO j=1,maxrp
      IF(sw.EQ.1)EXIT
!---------------------3.1------------------------!
!     Calculate traveltime gradient vector for current cell using
!     a first-order or second-order scheme.
!
      IF(igref.EQ.1)THEN
!
!        In this case, we are in the refined grid.
!
!        First order scheme applied here.
!
         dtx=ttnr(ipzr,ipxr+1)-ttnr(ipzr,ipxr)
         dtx=dtx+ttnr(ipzr+1,ipxr+1)-ttnr(ipzr+1,ipxr)
         dtx=dtx/(2.0*earth*dnxr)
         dtz=ttnr(ipzr+1,ipxr)-ttnr(ipzr,ipxr)
         dtz=dtz+ttnr(ipzr+1,ipxr+1)-ttnr(ipzr,ipxr+1)
         dtz=dtz/(2.0*earth*SIN(rgx(j))*dnzr)
      ELSE
!
!        Here, we are in the coarse grid.
!
!        First order scheme applied here.
!
         dtx=ttn(ipz,ipx+1)-ttn(ipz,ipx)
         dtx=dtx+ttn(ipz+1,ipx+1)-ttn(ipz+1,ipx)
         dtx=dtx/(2.0*earth*dnx)
         dtz=ttn(ipz+1,ipx)-ttn(ipz,ipx)
         dtz=dtz+ttn(ipz+1,ipx+1)-ttn(ipz,ipx+1)
         dtz=dtz/(2.0*earth*SIN(rgx(j))*dnz)
      ENDIF
!-----------------------3.2--------------------!
!     Calculate the next ray path point
!
      rd1=SQRT(dtx**2+dtz**2)
      rgx(j+1)=rgx(j)-dpl*dtx/(earth*rd1)
      rgz(j+1)=rgz(j)-dpl*dtz/(earth*SIN(rgx(j))*rd1)
!-----------------------3.3--------------------!
!     Determine which cell the new ray endpoint
!     lies in.
!
      ipxo=ipx
      ipzo=ipz
      IF(asgr.EQ.1)THEN
!
!        Here, we test to see whether the ray endpoint lies
!        within a cell of the refined grid
!
         ipxr=INT((rgx(j+1)-goxr)/dnxr)+1
         ipzr=INT((rgz(j+1)-gozr)/dnzr)+1
         igref=1
         IF(ipxr.LT.1.OR.ipxr.GE.nnxr)igref=0
         IF(ipzr.LT.1.OR.ipzr.GE.nnzr)igref=0
         IF(igref.EQ.1)THEN
            IF(nstsr(ipzr,ipxr).NE.0.OR.nstsr(ipzr+1,ipxr).NE.0)igref=0
            IF(nstsr(ipzr,ipxr+1).NE.0.OR.nstsr(ipzr+1,ipxr+1).NE.0)igref=0
         ENDIF
         ipx=INT((rgx(j+1)-gox)/dnx)+1
         ipz=INT((rgz(j+1)-goz)/dnz)+1
      ELSE
         ipx=INT((rgx(j+1)-gox)/dnx)+1
         ipz=INT((rgz(j+1)-goz)/dnz)+1
         igref=0
      ENDIF
!
!     Test the proximity of the source to the ray end point.
!     If it is less than dpl then we are done
!
      sred=((scx-rgx(j+1))*earth)**2
      sred=sred+((scz-rgz(j+1))*earth*SIN(rgx(j+1)))**2
      sred=SQRT(sred)
      sw=0
      IF(sred.LT.2.0*dpl)THEN
         rgx(j+2)=scx
         rgz(j+2)=scz
         nrp=j+2
         sw=1
!fang!         IF(cfd.NE.1)EXIT
      ENDIF
!
!     Due to the method for calculating traveltime gradient, if the
!     the ray end point lies in the source cell, then we are also done.
!
      IF(sw.EQ.0)THEN
         IF(asgr.EQ.1)THEN
            IF(igref.EQ.1)THEN
               IF(ipxr.EQ.isx)THEN
                  IF(ipzr.EQ.isz)THEN
                     rgx(j+2)=scx
                     rgz(j+2)=scz
                     nrp=j+2
                     sw=1
 !fang!                    IF(cfd.NE.1)EXIT
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF(ipx.EQ.isx)THEN
               IF(ipz.EQ.isz)THEN
                  rgx(j+2)=scx
                  rgz(j+2)=scz
                  nrp=j+2
                  sw=1
 !fang!                 IF(cfd.NE.1)EXIT
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     Test whether ray path segment extends beyond
!     box boundaries
!
      IF(ipx.LT.1)THEN
         rgx(j+1)=gox
         ipx=1
         rbint=1
      ENDIF
      IF(ipx.GE.nnx)THEN
         rgx(j+1)=gox+(nnx-1)*dnx
         ipx=nnx-1
         rbint=1
      ENDIF
      IF(ipz.LT.1)THEN
         rgz(j+1)=goz
         ipz=1
         rbint=1
      ENDIF
      IF(ipz.GE.nnz)THEN
         rgz(j+1)=goz+(nnz-1)*dnz
         ipz=nnz-1
         rbint=1
      ENDIF

! -----------------------------3.4---------------------------------!
!     Calculate the Frechet derivatives if required.
!--------------------------------------------------------------------!
!-------------------3.4.0--------------------------!
! Liu ! Calculate angle psi from [rgx(j),rgz(j)] to [rgx(j+1),rgz(j+1)]
 ! Must assume rgx is in colatitude rad
 !  azdist use latitude in degree, all variables are in degree
 ! fortran cos use rad
 ! error: azdist from event to station, so put the first point in the back.
 ! azdist(stalat, stalon, evtlat, evtlon, delta, az, baz)
 ! firt point is from reveiver.
   rgx1=(pi/2-rgx(j))*180.0/pi
   rgz1=rgz(j)*180.0/pi
   rgx2=(pi/2-rgx(j+1))*180.0/pi
   rgz2=rgz(j+1)*180.0/pi
   call azdist(rgx2,rgz2,rgx1,rgz1,delta,az,baz)
   rgaz(j)=az

   rgpsi=az/180*pi
  ! write(*,*)j,az,rgaz(j),rgpsi

 !fang!     IF(cfd.EQ.1)THEN
!-------------------3.4.1-------------------------!
!        First determine which B-spline cell the refined cells
!        containing the ray path segment lies in. If they lie
!        in more than one, then we need to divide the problem
!        into separate parts (up to three).
!
!       covert propagation grid coordinate to inversion grid.
         ivx=INT((ipx-1)/gdx)+1
         ivz=INT((ipz-1)/gdz)+1
         ivxo=INT((ipxo-1)/gdx)+1
         ivzo=INT((ipzo-1)/gdz)+1
!
!        Calculate up to two hit points between straight
!        ray segment and cell faces.
!
         nhp=0
         IF(ivx.NE.ivxo)THEN
            nhp=nhp+1
            IF(ivx.GT.ivxo)THEN
               xi=gox+(ivx-1)*dvx
            ELSE
               xi=gox+ivx*dvx
            ENDIF
            vrat(nhp)=(xi-rgx(j))/(rgx(j+1)-rgx(j))
            chp(nhp)=1
         ENDIF
         IF(ivz.NE.ivzo)THEN
            nhp=nhp+1
            IF(ivz.GT.ivzo)THEN
               zi=goz+(ivz-1)*dvz
            ELSE
               zi=goz+ivz*dvz
            ENDIF
            rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
            IF(nhp.EQ.1)THEN
               vrat(nhp)=rd1
               chp(nhp)=2
            ELSE
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=2
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=2
               ENDIF
            ENDIF
         ENDIF
         nhp=nhp+1
         vrat(nhp)=1.0
         chp(nhp)=0
!-------------3.4.2----------------!
!        Calculate the velocity, v and w values of the
!        first point [rgx(j),rgz(j)]
!
         drx=(rgx(j)-gox)-(ipxo-1)*dnx
         drz=(rgz(j)-goz)-(ipzo-1)*dnz
         vel=0.0
         DO l=1,2
            DO m=1,2
               produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
               produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
               IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx)THEN
                  vel=vel+veln(ipzo-1+m,ipxo-1+l)*produ
               ENDIF
            ENDDO
         ENDDO
         drx=(rgx(j)-gox)-(ivxo-1)*dvx
         drz=(rgz(j)-goz)-(ivzo-1)*dvz
         v=drx/dvx
         w=drz/dvz
!
!        Calculate the 12 basis values at the point
!
         vi(1)=(1.0-v)**3/6.0
         vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
         vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
         vi(4)=v**3/6.0
         wi(1)=(1.0-w)**3/6.0
         wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
         wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
         wi(4)=w**3/6.0
         ivxt=ivxo
         ivzt=ivzo
!-------------3.4.3-------------------!
!        Now loop through the one or more sub-segments of the
!        ray path segment and calculate partial derivatives
!
         DO k=1,nhp
            velo=vel
            vio=vi
            wio=wi
            IF(k.GT.1)THEN
               IF(chp(k-1).EQ.1)THEN
                  ivxt=ivx
               ELSE IF(chp(k-1).EQ.2)THEN
                  ivzt=ivz
               ENDIF
            ENDIF
!
!           Calculate the velocity, v and w values of the
!           new point
!
            rigz=rgz(j)+vrat(k)*(rgz(j+1)-rgz(j))
            rigx=rgx(j)+vrat(k)*(rgx(j+1)-rgx(j))
            ipxt=INT((rigx-gox)/dnx)+1
            ipzt=INT((rigz-goz)/dnz)+1
            drx=(rigx-gox)-(ipxt-1)*dnx
            drz=(rigz-goz)-(ipzt-1)*dnz
            vel=0.0
            DO m=1,2
               DO n=1,2
                  produ=(1.0-ABS(((n-1)*dnz-drz)/dnz))
                  produ=produ*(1.0-ABS(((m-1)*dnx-drx)/dnx))
                  IF(ipzt-1+n.LE.nnz.AND.ipxt-1+m.LE.nnx)THEN
                     vel=vel+veln(ipzt-1+n,ipxt-1+m)*produ
                  ENDIF
               ENDDO
            ENDDO
            drx=(rigx-gox)-(ivxt-1)*dvx
            drz=(rigz-goz)-(ivzt-1)*dvz
            v=drx/dvx
            w=drz/dvz
!
!           Calculate the 8 basis values at the new point
!
            vi(1)=(1.0-v)**3/6.0
            vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
            vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
            vi(4)=v**3/6.0
            wi(1)=(1.0-w)**3/6.0
            wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
            wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
            wi(4)=w**3/6.0
!
!           Calculate the incremental path length
!
            IF(k.EQ.1)THEN
               dinc=vrat(k)*dpl
            ELSE
               dinc=(vrat(k)-vrat(k-1))*dpl
            ENDIF
!
!           Now compute the 16 contributions to the partial
!           derivatives.
!
            DO l=1,4
               DO m=1,4
 ! isotropy part
                  rdc1=vi(m)*wi(l)/vel**2
                  rdc2=vio(m)*wio(l)/velo**2
                  rd1= -(rdc1+rdc2)*dinc/2.0
                  rd2=fdm(ivzt-2+l,ivxt-2+m)
                  fdm(ivzt-2+l,ivxt-2+m)=rd1+rd2
! anisotropy part cos--Liu----error fdm
                  rd1= -(rdc1*cos(2.0*rgpsi)+rdc2*cos(2.0*rgpsi))*dinc/2.0
                  rd2=fdmc(ivzt-2+l,ivxt-2+m)
                  fdmc(ivzt-2+l,ivxt-2+m)=rd1+rd2
! anisotropy part sin
                  rd1= -(rdc1*sin(2.0*rgpsi)+rdc2*sin(2.0*rgpsi))*dinc/2.0
                  rd2=fdms(ivzt-2+l,ivxt-2+m)
                  fdms(ivzt-2+l,ivxt-2+m)=rd1+rd2
               ENDDO
            ENDDO
         ENDDO
 !fang!     ENDIF
!fang!      IF(j.EQ.maxrp.AND.sw.EQ.0)THEN
!fang!         WRITE(6,*)'Error with ray path detected!!!'
!fang!         WRITE(6,*)'Source id: ',csid
!fang!         WRITE(6,*)'Receiver id: ',i
!fang!      ENDIF
   ENDDO
!
!  Write ray paths to output file
!
!fang!   IF(wrgf.EQ.csid.OR.wrgf.LT.0)THEN
  IF(writepath) THEN
!      WRITE(40,*)'>',nrp
      WRITE(40,'(a,f4.1)')'>',Tperiod
      DO j=1,nrp
         rayx=(pi/2-rgx(j))*180.0/pi
         rayz=rgz(j)*180.0/pi
         WRITE(40,*)rayz,rayx
      ENDDO
  ENDIF

!  open(88,file='Ray_Angle.dat',position='Append',action='write')
!      WRITE(88,'(a,f4.1)')'>',Tperiod
!      DO j=1,nrp
!         rayx=(pi/2-rgx(j))*180.0/pi
!         rayz=rgz(j)*180.0/pi
!         WRITE(88,*)j,rayz,rayx,rgaz(j)
!      ENDDO
!      close(88)

!     open(88,file='Ray_fdmc.dat',position='Append',action='write')
!     WRITE(88,'(a,f4.1)')'>',Tperiod
!     write(88,*)" jj    kk           fdm           fdmc           fdms        angle      cos(2*psi)       sin(2*psi)"
!     do jj=1,nvz
!           do kk=1,nvx
!                 if(abs(fdm(jj,kk)).ge.ftol) then
!                    write(88,*)jj,kk,fdm(jj,kk),fdmc(jj,kk),fdms(jj,kk),sum(rgaz(1:nrp))/nrp,&
!                    cos(2*sum(rgaz(1:nrp))/nrp/180*pi),sin(2*sum(rgaz(1:nrp))/nrp/180*pi)
!                endif
!            enddo
!       enddo
!     close(88)

!fang!   ENDIF
!
!  Write partial derivatives to output file
!
!fang!   IF(cfd.EQ.1)THEN
!fang!!
!fang!!     Determine the number of non-zero elements.
!fang!!
!fang!      isum=0
!fang!      DO j=0,nvz+1
!fang!         DO k=0,nvx+1
!fang!            IF(ABS(fdm(j,k)).GE.ftol)isum=isum+1
!fang!         ENDDO
!fang!      ENDDO
!fang!      WRITE(50)isum
!fang!      isum=0
!fang!      DO j=0,nvz+1
!fang!         DO k=0,nvx+1
!fang!            isum=isum+1
!fang!            IF(ABS(fdm(j,k)).GE.ftol)WRITE(50)isum,fdm(j,k)
!fang!         ENDDO
!fang!      ENDDO
!fang!   ENDIF
!fang!ENDDO
!fang!IF(cfd.EQ.1)THEN
!fang!   DEALLOCATE(fdm, STAT=checkstat)
!fang!   IF(checkstat > 0)THEN
!fang!      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: fdm'
!fang!   ENDIF
!fang!ENDIF
DEALLOCATE(rgx,rgz, STAT=checkstat)
DEALLOCATE(rgaz, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpathsAzim: rgx,rgz'
ENDIF
END SUBROUTINE rpathsAzim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE azdist(stalat, stalon, evtlat, evtlon,&
                    delta, az, baz)
!
! Subroutine to calculate the Great Circle Arc distance
!    between two sets of geographic coordinates
!
! Given:  stalat => Latitude of first point (+N, -S) in degrees
!	  stalon => Longitude of first point (+E, -W) in degrees
!	  evtlat => Latitude of second point
!	  evtlon => Longitude of second point
!
! Returns:  delta => Great Circle Arc distance in degrees
!	    az    => Azimuth from pt. 1 to pt. 2 in degrees
!	    baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
!
! If you are calculating station-epicenter pairs, pt. 1 is the station
!
! Equations take from Bullen, pages 154, 155
!
! T. Owens, September 19, 1991
!           Sept. 25 -- fixed az and baz calculations
!           Dec. 2006, changed for fortran95
!           May, 2007 -- added predel to get around OSX acos round-off NaN issue
!
      double precision scolat, slon, ecolat, elon
      double precision a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk
      double precision rhs1,rhs2,sph,rad,del,daz,dbaz,pi
!
      pi=3.1415926535898
      piby2=pi/2.
      rad=2.*pi/360.
!
! scolat and ecolat are the geocentric colatitudes
! as defined by Richter (pg. 318)
!
! Earth Flattening of 1/298.257 take from Bott (pg. 3)
!
      sph=1.0/298.257
!
      scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(stalat)*rad))
      ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(evtlat)*rad))
      slon=dble(stalon)*rad
      elon=dble(evtlon)*rad
!
!  a - e are as defined by Bullen (pg. 154, Sec 10.2)
!     These are defined for the pt. 1
!
      a=sin(scolat)*cos(slon)
      b=sin(scolat)*sin(slon)
      c=cos(scolat)
      d=sin(slon)
      e=-cos(slon)
      g=-c*e
      h=c*d
      k=-sin(scolat)
!
!  aa - ee are the same as a - e, except for pt. 2
!
      aa=sin(ecolat)*cos(elon)
      bb=sin(ecolat)*sin(elon)
      cc=cos(ecolat)
      dd=sin(elon)
      ee=-cos(elon)
      gg=-cc*ee
      hh=cc*dd
      kk=-sin(ecolat)
!
!  Bullen, Sec 10.2, eqn. 4
!
      predel=a*aa + b*bb + c*cc
      if(abs(predel+1.).lt..000001) then
        predel=-1.
      endif
      if(abs(predel-1.).lt..000001) then
        predel=1.
      endif
      del=acos(predel)
      delta=del/rad
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 1 is unprimed, so this is technically the baz
!
!  Calculate baz this way to avoid quadrant problems
!
      rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.
      rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.
      dbaz=atan2(rhs1,rhs2)
      if(dbaz.lt.0.0d0) dbaz=dbaz+2*pi
      baz=dbaz/rad
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 2 is unprimed, so this is technically the az
!
      rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.
      rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.
      daz=atan2(rhs1,rhs2)
      if(daz.lt.0.0d0) daz=daz+2*pi
      az=daz/rad
!
!   Make sure 0.0 is always 0.0, not 360.
!
      if(abs(baz-360.).lt..00001) baz=0.0
      if(abs(az-360.).lt..00001) az=0.0
      return
END SUBROUTINE azdist


