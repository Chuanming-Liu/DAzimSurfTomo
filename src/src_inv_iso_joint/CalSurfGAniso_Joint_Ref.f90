!---------------------------------------------------------------------------------!
! package for the first iteration-which not use the Vs+dv
!---------------------------------------------------------------------------------!
      subroutine CalRayleighPhase(nx,ny,nz,vel,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
        use omp_lib
        implicit none
! INPUT
        integer nx,ny,nz
        real vel(nx,ny,nz)

        integer iwave,igr
        real minthk
        real depz(nz)
        integer kmaxRc
        real*8 tRc(kmaxRc)
! OUTPUT
        real*8 pvRc(nx*ny,kmaxRc)
! PARAMETETER LIST

        real vpz(nz),vsz(nz),rhoz(nz)
	integer mmax,iflsph,mode,rmax
        integer ii,jj,k,i,nn,kk
    	integer,parameter::NL=200
    	integer,parameter::NP=60
        real*8 cg1(NP),cg2(NP),cga,cgRc(NP)
    	real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    	real depm(NL),vpm(NL),vsm(NL),rhom(NL),thkm(NL)
	real dlnVs,dlnVp,dlnrho

        mmax=nz
        iflsph=1
        mode=1
       pvRc=0
	!print*,'depth kernel begin...'
!$omp parallel &
!$omp default(private) &
!$omp shared(depz,nx,ny,nz,minthk,kmaxRc,mmax,vel) &
!$omp shared(tRc,pvRc,iflsph,iwave,mode,igr)
!$omp do
        do jj=1,ny
    	do ii=1,nx
    	vsz(1:nz)=vel(ii,jj,1:nz)
        ! some other emperical relationship maybe better,
    	do k=1,nz
        vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
        0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
    	rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
               0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
               0.000106*vpz(k)**5
        enddo

    	call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax,rdep,&
        rvp,rvs,rrho,rthk)
    	call surfdisp96(rthk,rvp,rvs,rrho,rmax,iflsph,iwave,mode,igr,kmaxRc,&
        tRc,cgRc)
        pvRc((jj-1)*nx+ii,1:kmaxRc)=cgRc(1:kmaxRc)
        !print*,cgRc(1:kmaxRc)
    	enddo
    	enddo
!$omp end do
!$omp end parallel

   end subroutine CalRayleighPhase



subroutine refineGrid2LayerMdl(minthk0,mmax,dep,vp,vs,rho,&
                  rmax,rdep,rvp,rvs,rrho,rthk)
!!--------------------------------------------------------------------c
!!refine grid based model to layerd based model
!!INPUT:   minthk: is the minimum thickness of the refined layered model
!!         mmax: number of depth grid points in the model
!!         dep, vp, vs, rho: the depth-grid model parameters
!!         rmax: number of layers in the fined layered model
!!         rdep, rvp, rvs, rrho, rthk: the refined layered velocity model
!!
        implicit none
        integer NL
        parameter (NL=200)
        integer mmax,rmax
        real minthk0
        real minthk
        real dep(*),vp(*),vs(*),rho(*)
        real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
        integer nsublay(NL)
        real thk,newthk,initdep
        integer i,j,k,ngrid

        k = 0
        initdep = 0.0
        do i = 1, mmax-1
           thk = dep(i+1)-dep(i)
	   minthk = thk/minthk0
           nsublay(i) = int((thk+1.0e-4)/minthk) + 1
           ngrid = nsublay(i)+1
           newthk = thk/nsublay(i)
           do j = 1, nsublay(i)
              k = k + 1
              rthk(k) = newthk
              rdep(k) = initdep + rthk(k)
              initdep = rdep(k)
              rvp(k) = vp(i)+(2*j-1)*(vp(i+1)-vp(i))/(2*nsublay(i))
              rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*nsublay(i))
              rrho(k) = rho(i)+(2*j-1)*(rho(i+1)-rho(i))/(2*nsublay(i))
           enddo
        enddo
!! half space model
        k = k + 1
        rthk(k) = 0.0
        rvp(k) = vp(mmax)
        rvs(k) = vs(mmax)
        rrho(k) = rho(mmax)
	rdep(k) = dep(mmax)

        rmax = k

!!       do i = 1, mmax
!!          write(*,*) dep(i),vp(i),vs(i),rho(i)
!!       enddo
!!       print *, '---------------------------------'
!!       do i = 1, rmax
!!          write(*,*) rdep(i),rthk(i),rvp(i),rvs(i),rrho(i)
!!       enddo

        return
        end subroutine

!---------------------------------------------------------------------------------!
! package for the  iteration-which use the Vs+dv
!---------------------------------------------------------------------------------!

subroutine CalRayPhIteration(nx,ny,nz,vel,vsRela,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
        use omp_lib
        implicit none
! INPUT
        integer nx,ny,nz
        real vel(nx,ny,nz)
        REAL vsRela(nx-2,ny-2,nz-1)

        integer iwave,igr
        real minthk
        real depz(nz)
        integer kmaxRc
        real*8 tRc(kmaxRc)
! OUTPUT
        real*8 pvRc(nx*ny,kmaxRc)
! PARAMETETER LIST

        real vpz(nz),vsz(nz),rhoz(nz)
	integer mmax,iflsph,mode,rmax
        integer ii,jj,k,i,nn,kk
                INTEGER nvx,nvz
    	integer,parameter::NL=200
    	integer,parameter::NP=60
        real*8 cg1(NP),cg2(NP),cga,cgRc(NP)
!    	real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    	real depm(NL),vpm(NL),vsm(NL),rhom(NL),thkm(NL)
	real dlnVs,dlnVp,dlnrho
! PARAMETER ITERATION
        REAL dvs(nx,ny,nz-1)
        REAL vsPt(nz-1)
        REAL rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
        REAl rvsPt(NL), rvpPt(NL),rhoPt(NL)
! handle outside point
        dvs=0
        dvs(2:nx-1,2:ny-1,1:nz-1)=vsRela(1:nx-2,1:ny-2,1:nz-1)

        nvx=nx-2
        nvz=ny-2
        iwave=2

        iflsph=1
        mode=1
        pvRc=0
	!print*,'depth kernel begin...'
!$omp parallel &
!$omp default(private) &
!$omp shared(depz,nx,ny,nz,minthk,kmaxRc,vel,dvs) &
!$omp shared(tRc,pvRc,iflsph,iwave,mode,igr)
!$omp do
        do jj=1,ny
    	do ii=1,nx
               vsPt(1:nz-1)=dvs(ii,jj,1:nz-1)
    	       vsz(1:nz)=vel(ii,jj,1:nz)

    	       call refineGrid2LayerIte(minthk,nz,depz,vsz,vsPt,rmax,rdep,rvsPt,rthk)

             do k=1,rmax
                    rvpPt(k)=0.9409 + 2.0947*rvsPt(k) - 0.8206*rvsPt(k)**2+ &
                             0.2683*rvsPt(k)**3 - 0.0251*rvsPt(k)**4
    	            rhoPt(k)=1.6612*rvpPt(k) - 0.4721*rvpPt(k)**2 + &
                        0.0671*rvpPt(k)**3 - 0.0043*rvpPt(k)**4 + &
                       0.000106*rvpPt(k)**5
             enddo

            call surfdisp96(rthk,rvpPt,rvsPt,rhoPt,rmax,iflsph,iwave,mode,igr,kmaxRc,&
            tRc,cgRc)
           if(minval(cgRc(1:kmaxRc)).LT.0.0) STOP "Fwd Phase Error!"

           pvRc((jj-1)*nx+ii,1:kmaxRc)=cgRc(1:kmaxRc)
        enddo
    	enddo
!$omp end do
!$omp end parallel




end subroutine


subroutine refineGrid2LayerIte(minthk0,mmax,dep,vs,vsPt,rmax,rdep,rvsPt,rthk)
       IMPLICIT NONE
        integer NL
        parameter (NL=200)
! INPUT
        REAL minthk0
        INTEGER mmax
        REAL dep(mmax),vs(mmax),vsPt(mmax-1)
! OUTPUT
        REAL rdep(NL),rvs(NL),rthk(NL),rvsPt(NL)
        INTEGER rmax
!  PARAMETER
        real minthk
        integer nsublay(NL)
        real thk,newthk,initdep
        integer i,j,k,ngrid
!--------------------------------------------------------------------!
! refine grid based model to layerd based model
!  minthk0: number of sub layer
!  mmax: number of depth grid points in the model---->nz
!  dep, vp, vs, rho: the depth-grid model parameters
!  rmax: number of layers in the fined layered model
!  rdep, rvp, rvs, rrho, rthk: the refined layered velocity model
! rvs: refined grid vs model
! rvsPt: rvs+ inversion result, which is same at i th layer.
!---------------------------------------------------------------------!
        k = 0
        initdep = 0.0
        rvs=0
        rvsPt=0
        do i = 1, mmax-1
           thk = dep(i+1)-dep(i)
	   minthk = thk/minthk0
           nsublay(i) = int((thk+1.0e-4)/minthk) + 1
           ngrid = nsublay(i)+1
           newthk = thk/nsublay(i)
           do j = 1, nsublay(i)
              k = k + 1
              rthk(k) = newthk
              rdep(k) = initdep + rthk(k)
              initdep = rdep(k)
              rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*nsublay(i))
              rvsPt(k)=rvs(k)+rvs(k)*vsPt(i)
           enddo
        enddo
!! half space model
        k = k + 1
        rthk(k) = 0.0
        rvs(k) = vs(mmax)
        rvsPt(k)=vs(mmax)
	rdep(k) = dep(mmax)
        rmax = k

        return
        end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(6)
INTEGER :: checkstat
INTEGER, SAVE :: nvx,nvz,nnx,nnz,fom,gdx,gdz
INTEGER, SAVE :: vnl,vnr,vnt,vnb,nrnx,nrnz,sgdl,rbint
INTEGER, SAVE :: nnxr,nnzr,asgr
INTEGER, DIMENSION (:,:), ALLOCATABLE :: nsts,nstsr,srs
REAL(KIND=i10), SAVE :: gox,goz,dnx,dnz,dvx,dvz,snb,earth
REAL(KIND=i10), SAVE :: goxd,gozd,dvxd,dvzd,dnxd,dnzd
REAL(KIND=i10), SAVE :: drnx,drnz,gorx,gorz
REAL(KIND=i10), SAVE :: dnxr,dnzr,goxr,gozr
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: velv,veln,velnb
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttn,ttnr
!REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: rcx,rcz
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
!!!--------------------------------------------------------------
!!	modified by Hongjian Fang @ USTC
!	real,dimension(:),allocatable,save::rw
!	integer,dimension(:),allocatable,save::iw,col
!	real,dimension(:,:,:),allocatable::vpf,vsf
!	real,dimension(:),allocatable,save::obst,cbst,wt,dtres
!!	integer,dimension(:),allocatable,save::cbst_stat
!	real,dimension(:,:,:),allocatable,save::sen_vs,sen_vp,sen_rho
!!!	real,dimension(:,:,:),allocatable,save::sen_vsRc,sen_vpRc,sen_rhoRc
!!!	real,dimension(:,:,:),allocatable,save::sen_vsRg,sen_vpRg,sen_rhoRg
!!!	real,dimension(:,:,:),allocatable,save::sen_vsLc,sen_vpLc,sen_rhoLc
!!!	real,dimension(:,:,:),allocatable,save::sen_vsLg,sen_vpLg,sen_rhoLg
!!!	integer,save:: count1,count2
!	integer*8,save:: nar
!	integer,save:: iter,maxiter
!!!--------------------------------------------------------------
!
! nvx,nvz = B-spline vertex values
! dvx,dvz = B-spline vertex separation
! velv(i,j) = velocity values at control points
! nnx,nnz = Number of nodes of grid in x and z
! nnxr,nnzr = Number of nodes of refined grid in x and z
! gox,goz = Origin of grid (theta,phi)
! goxr, gozr = Origin of refined grid (theta,phi)
! dnx,dnz = Node separation of grid in  x and z
! dnxr,dnzr = Node separation of refined grid in x and z
! veln(i,j) = velocity values on a refined grid of nodes
! velnb(i,j) = Backup of veln required for source grid refinement
! ttn(i,j) = traveltime field on the refined grid of nodes
! ttnr(i,j) = ttn for refined grid
! nsts(i,j) = node status (-1=far,0=alive,>0=close)
! nstsr(i,j) = nsts for refined grid
! checkstat = check status of memory allocation
! fom = use first-order(0) or mixed-order(1) scheme
! snb = Maximum size of narrow band as fraction of nnx*nnz
! nrc = number of receivers
! rcx(i),rcz(i) = (x,z) coordinates of receivers
! earth = radius of Earth (in km)
! goxd,gozd = gox,goz in degrees
! dvxd,dvzd = dvx,dvz in degrees
! dnzd,dnzd = dnx,dnz in degrees
! gdx,gdz = grid dicing in x and z
! vnl,vnr,vnb,vnt = Bounds of refined grid
! nrnx,nrnz = Number of nodes in x and z for refined grid
! gorx,gorz = Grid origin of refined grid
! sgdl = Source grid dicing level
! rbint = Ray-boundary intersection (0=no, 1=yes).
! asgr = Apply source grid refinement (0=no,1=yes)
! srs = Source-receiver status (0=no path, 1=path exists)
!
END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module contains all the subroutines used to calculate
! the first-arrival traveltime field through the grid.
! Subroutines are:
! (1) travel
! (2) fouds1
! (3) fouds2
! (4) addtree
! (5) downtree
! (6) updtree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE traveltime
USE globalp
IMPLICIT NONE
INTEGER ntr
TYPE backpointer
   INTEGER(KIND=2) :: px,pz
END TYPE backpointer
TYPE(backpointer), DIMENSION (:), ALLOCATABLE :: btg
!
! btg = backpointer to relate grid nodes to binary tree entries
! px = grid-point in x
! pz = grid-point in z
! ntr = number of entries in binary tree
!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE travel(scx,scz,urg)
IMPLICIT NONE
INTEGER :: isx,isz,sw,i,j,ix,iz,urg,swrg
REAL(KIND=i10) :: scx,scz,vsrc,dsx,dsz,ds
REAL(KIND=i10), DIMENSION (2,2) :: vss
! isx,isz = grid cell indices (i,j,k) which contains source
! scx,scz = (r,x,y) location of source
! sw = a switch (0=off,1=on)
! ix,iz = j,k position of "close" point with minimum traveltime
! maxbt = maximum size of narrow band binary tree
! rd2,rd3 = substitution variables
! vsrc = velocity at source
! vss = velocity at nodes surrounding source
! dsx, dsz = distance from source to cell boundary in x and z
! ds = distance from source to nearby node
! urg = use refined grid (0=no,1=yes,2=previously used)
! swrg = switch to end refined source grid computation
!
! The first step is to find out where the source resides
! in the grid of nodes. The cell in which it resides is
! identified by the "north-west" node of the cell. If the
! source lies on the edge or corner (a node) of the cell, then
! this scheme still applies.
!
isx=INT((scx-gox)/dnx)+1
isz=INT((scz-goz)/dnz)+1
sw=0
IF(isx.lt.1.or.isx.gt.nnx)sw=1
IF(isz.lt.1.or.isz.gt.nnz)sw=1
IF(sw.eq.1)then
   scx=90.0-scx*180.0/pi
   scz=scz*180.0/pi
   WRITE(6,*)"Source lies outside bounds of model (lat,long)= ",scx,scz
   WRITE(6,*)"TERMINATING PROGRAM!!!"
   STOP
ENDIF
IF(isx.eq.nnx)isx=isx-1
IF(isz.eq.nnz)isz=isz-1
!
! Set all values of nsts to -1 if beginning from a source
! point.
!
IF(urg.NE.2)nsts=-1
!
! set initial size of binary tree to zero
!
ntr=0
IF(urg.EQ.2)THEN
!
!  In this case, source grid refinement has been applied, so
!  the initial narrow band will come from resampling the
!  refined grid.
!
   DO i=1,nnx
      DO j=1,nnz
         IF(nsts(j,i).GT.0)THEN
            CALL addtree(j,i)
         ENDIF
      ENDDO
   ENDDO
ELSE
!
!  In general, the source point need not lie on a grid point.
!  Bi-linear interpolation is used to find velocity at the
!  source point.
!
   nsts=-1
   DO i=1,2
      DO j=1,2
         vss(i,j)=veln(isz-1+j,isx-1+i)
      ENDDO
   ENDDO
   dsx=(scx-gox)-(isx-1)*dnx
   dsz=(scz-goz)-(isz-1)*dnz
   CALL bilinear(vss,dsx,dsz,vsrc)
!
!  Now find the traveltime at the four surrounding grid points. This
!  is calculated approximately by assuming the traveltime from the
!  source point to each node is equal to the the distance between
!  the two points divided by the average velocity of the points
!
   DO i=1,2
      DO j=1,2
         ds=SQRT((dsx-(i-1)*dnx)**2+(dsz-(j-1)*dnz)**2)
         ttn(isz-1+j,isx-1+i)=2.0*ds/(vss(i,j)+vsrc)
         CALL addtree(isz-1+j,isx-1+i)
      ENDDO
   ENDDO
ENDIF
!
! Now calculate the first-arrival traveltimes at the
! remaining grid points. This is done via a loop which
! repeats the procedure of finding the first-arrival
! of all "close" points, adding it to the set of "alive"
! points and updating the points surrounding the new "alive"
! point. The process ceases when the binary tree is empty,
! in which case all grid points are "alive".
!
DO WHILE(ntr.gt.0)
!
! First, check whether source grid refinement is
! being applied; if so, then there is a special
! exit condition.
!
IF(urg.EQ.1)THEN
   ix=btg(1)%px
   iz=btg(1)%pz
   swrg=0
   IF(ix.EQ.1)THEN
      IF(vnl.NE.1)swrg=1
   ENDIF
   IF(ix.EQ.nnx)THEN
      IF(vnr.NE.nnx)swrg=1
   ENDIF
   IF(iz.EQ.1)THEN
      IF(vnt.NE.1)swrg=1
   ENDIF
   IF(iz.EQ.nnz)THEN
      IF(vnb.NE.nnz)swrg=1
   ENDIF
   IF(swrg.EQ.1)THEN
      nsts(iz,ix)=0
      EXIT
   ENDIF
ENDIF
!
! Set the "close" point with minimum traveltime
! to "alive"
!
   ix=btg(1)%px
   iz=btg(1)%pz
   nsts(iz,ix)=0
!
! Update the binary tree by removing the root and
! sweeping down the tree.
!
   CALL downtree
!
! Now update or find values of up to four grid points
! that surround the new "alive" point.
!
! Test points that vary in x
!
   DO i=ix-1,ix+1,2
      IF(i.ge.1.and.i.le.nnx)THEN
         IF(nsts(iz,i).eq.-1)THEN
!
! This option occurs when a far point is added to the list
! of "close" points
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,i)
            ELSE
               CALL fouds2(iz,i)
            ENDIF
            CALL addtree(iz,i)
         ELSE IF(nsts(iz,i).gt.0)THEN
!
! This happens when a "close" point is updated
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,i)
            ELSE
               CALL fouds2(iz,i)
            ENDIF
            CALL updtree(iz,i)
         ENDIF
      ENDIF
   ENDDO
!
! Test points that vary in z
!
   DO i=iz-1,iz+1,2
      IF(i.ge.1.and.i.le.nnz)THEN
         IF(nsts(i,ix).eq.-1)THEN
!
! This option occurs when a far point is added to the list
! of "close" points
!
            IF(fom.eq.0)THEN
               CALL fouds1(i,ix)
            ELSE
               CALL fouds2(i,ix)
            ENDIF
            CALL addtree(i,ix)
         ELSE IF(nsts(i,ix).gt.0)THEN
!
! This happens when a "close" point is updated
!
            IF(fom.eq.0)THEN
               CALL fouds1(i,ix)
            ELSE
               CALL fouds2(i,ix)
            ENDIF
            CALL updtree(i,ix)
         ENDIF
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE travel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates a trial first-arrival traveltime
! at a given node from surrounding nodes using the
! First-Order Upwind Difference Scheme (FOUDS) of
! Sethian and Popovici (1999).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fouds1(iz,ix)
IMPLICIT NONE
INTEGER :: j,k,ix,iz,tsw1,swsol
REAL(KIND=i10) :: trav,travm,slown,tdsh,tref
REAL(KIND=i10) :: a,b,c,u,v,em,ri,risti
REAL(KIND=i10) :: rd1
!
! ix = NS position of node coordinate for determination
! iz = EW vertical position of node coordinate for determination
! trav = traveltime calculated for trial node
! travm = minimum traveltime calculated for trial node
! slown = slowness at (iz,ix)
! tsw1 = traveltime switch (0=first time,1=previously)
! a,b,c,u,v,em = Convenience variables for solving quadratic
! tdsh = local traveltime from neighbouring node
! tref = reference traveltime at neighbouring node
! ri = Radial distance
! risti = ri*sin(theta) at point (iz,ix)
! rd1 = dummy variable
! swsol = switch for solution (0=no solution, 1=solution)
!
! Inspect each of the four quadrants for the minimum time
! solution.
!
tsw1=0
slown=1.0/veln(iz,ix)
ri=earth
risti=ri*sin(gox+(ix-1)*dnx)
DO j=ix-1,ix+1,2
   DO k=iz-1,iz+1,2
      IF(j.GE.1.AND.j.LE.nnx)THEN
         IF(k.GE.1.AND.k.LE.nnz)THEN
!
!           There are seven solution options in
!           each quadrant.
!
            swsol=0
            IF(nsts(iz,j).EQ.0)THEN
               swsol=1
               IF(nsts(k,ix).EQ.0)THEN
                  u=ri*dnx
                  v=risti*dnz
                  em=ttn(k,ix)-ttn(iz,j)
                  a=u**2+v**2
                  b=-2.0*u**2*em
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
               ELSE
                  a=1.0
                  b=0.0
                  c=-slown**2*ri**2*dnx**2
                  tref=ttn(iz,j)
               ENDIF
            ELSE IF(nsts(k,ix).EQ.0)THEN
               swsol=1
               a=1.0
               b=0.0
               c=-(slown*risti*dnz)**2
               tref=ttn(k,ix)
            ENDIF
!
!           Now find the solution of the quadratic equation
!
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=tref+tdsh
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
                ENDIF
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO
ttn(iz,ix)=travm
END SUBROUTINE fouds1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates a trial first-arrival traveltime
! at a given node from surrounding nodes using the
! Mixed-Order (2nd) Upwind Difference Scheme (FOUDS) of
! Popovici and Sethian (2002).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fouds2(iz,ix)
IMPLICIT NONE
INTEGER :: j,k,j2,k2,ix,iz,tsw1
INTEGER :: swj,swk,swsol
REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
REAL(KIND=i10) :: a,b,c,u,v,em,ri,risti,rd1
!
! ix = NS position of node coordinate for determination
! iz = EW vertical position of node coordinate for determination
! trav = traveltime calculated for trial node
! travm = minimum traveltime calculated for trial node
! slown = slowness at (iz,ix)
! tsw1 = traveltime switch (0=first time,1=previously)
! a,b,c,u,v,em = Convenience variables for solving quadratic
! tdsh = local traveltime from neighbouring node
! tref = reference traveltime at neighbouring node
! ri = Radial distance
! risti = ri*sin(theta) at point (iz,ix)
! swj,swk = switches for second order operators
! tdiv = term to divide tref by depending on operator order
! swsol = switch for solution (0=no solution, 1=solution)
!
! Inspect each of the four quadrants for the minimum time
! solution.
!
tsw1=0
slown=1.0/veln(iz,ix)
ri=earth
risti=ri*sin(gox+(ix-1)*dnx)
DO j=ix-1,ix+1,2
   IF(j.GE.1.AND.j.LE.nnx)THEN
      swj=-1
      IF(j.eq.ix-1)THEN
         j2=j-1
         IF(j2.GE.1)THEN
            IF(nsts(iz,j2).EQ.0)swj=0
         ENDIF
      ELSE
         j2=j+1
         IF(j2.LE.nnx)THEN
            IF(nsts(iz,j2).EQ.0)swj=0
         ENDIF
      ENDIF
      IF(nsts(iz,j).EQ.0.AND.swj.EQ.0)THEN
         swj=-1
         IF(ttn(iz,j).GT.ttn(iz,j2))THEN
            swj=0
         ENDIF
      ELSE
         swj=-1
      ENDIF
      DO k=iz-1,iz+1,2
         IF(k.GE.1.AND.k.LE.nnz)THEN
            swk=-1
            IF(k.eq.iz-1)THEN
               k2=k-1
               IF(k2.GE.1)THEN
                  IF(nsts(k2,ix).EQ.0)swk=0
               ENDIF
            ELSE
               k2=k+1
               IF(k2.LE.nnz)THEN
                  IF(nsts(k2,ix).EQ.0)swk=0
               ENDIF
            ENDIF
            IF(nsts(k,ix).EQ.0.AND.swk.EQ.0)THEN
               swk=-1
               IF(ttn(k,ix).GT.ttn(k2,ix))THEN
                  swk=0
               ENDIF
            ELSE
               swk=-1
            ENDIF
!
!           There are 8 solution options in
!           each quadrant.
!
            swsol=0
            IF(swj.EQ.0)THEN
               swsol=1
               IF(swk.EQ.0)THEN
                  u=2.0*ri*dnx
                  v=2.0*risti*dnz
                  em=4.0*ttn(iz,j)-ttn(iz,j2)-4.0*ttn(k,ix)
                  em=em+ttn(k2,ix)
                  a=v**2+u**2
                  b=2.0*em*u**2
                  c=u**2*(em**2-slown**2*v**2)
                  tref=4.0*ttn(iz,j)-ttn(iz,j2)
                  tdiv=3.0
               ELSE IF(nsts(k,ix).EQ.0)THEN
                  u=risti*dnz
                  v=2.0*ri*dnx
                  em=3.0*ttn(k,ix)-4.0*ttn(iz,j)+ttn(iz,j2)
                  a=v**2+9.0*u**2
                  b=6.0*em*u**2
                  c=u**2*(em**2-slown**2*v**2)
                  tref=ttn(k,ix)
                  tdiv=1.0
               ELSE
                  u=2.0*ri*dnx
                  a=1.0
                  b=0.0
                  c=-u**2*slown**2
                  tref=4.0*ttn(iz,j)-ttn(iz,j2)
                  tdiv=3.0
               ENDIF
            ELSE IF(nsts(iz,j).EQ.0)THEN
               swsol=1
               IF(swk.EQ.0)THEN
                  u=ri*dnx
                  v=2.0*risti*dnz
                  em=3.0*ttn(iz,j)-4.0*ttn(k,ix)+ttn(k2,ix)
                  a=v**2+9.0*u**2
                  b=6.0*em*u**2
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
                  tdiv=1.0
               ELSE IF(nsts(k,ix).EQ.0)THEN
                  u=ri*dnx
                  v=risti*dnz
                  em=ttn(k,ix)-ttn(iz,j)
                  a=u**2+v**2
                  b=-2.0*u**2*em
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
                  tdiv=1.0
               ELSE
                  a=1.0
                  b=0.0
                  c=-slown**2*ri**2*dnx**2
                  tref=ttn(iz,j)
                  tdiv=1.0
               ENDIF
            ELSE
               IF(swk.EQ.0)THEN
                  swsol=1
                  u=2.0*risti*dnz
                  a=1.0
                  b=0.0
                  c=-u**2*slown**2
                  tref=4.0*ttn(k,ix)-ttn(k2,ix)
                  tdiv=3.0
               ELSE IF(nsts(k,ix).EQ.0)THEN
                  swsol=1
                  a=1.0
                  b=0.0
                  c=-slown**2*risti**2*dnz**2
                  tref=ttn(k,ix)
                  tdiv=1.0
               ENDIF
            ENDIF
!
!           Now find the solution of the quadratic equation
!
            IF(swsol.EQ.1)THEN
               rd1=b**2-4.0*a*c
               IF(rd1.LT.0.0)rd1=0.0
               tdsh=(-b+sqrt(rd1))/(2.0*a)
               trav=(tref+tdsh)/tdiv
               IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
               ELSE
                  travm=trav
                  tsw1=1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF
ENDDO
ttn(iz,ix)=travm
END SUBROUTINE fouds2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine adds a value to the binary tree by
! placing a value at the bottom and pushing it up
! to its correct position.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE addtree(iz,ix)
IMPLICIT NONE
INTEGER :: ix,iz,tpp,tpc
TYPE(backpointer) :: exch
!
! ix,iz = grid position of new addition to tree
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
!
! First, increase the size of the tree by one.
!
ntr=ntr+1
!
! Put new value at base of tree
!
nsts(iz,ix)=ntr
btg(ntr)%px=ix
btg(ntr)%pz=iz
!
! Now filter the new value up to its correct position
!
tpc=ntr
tpp=tpc/2
DO WHILE(tpp.gt.0)
   IF(ttn(iz,ix).lt.ttn(btg(tpp)%pz,btg(tpp)%px))THEN
      nsts(iz,ix)=tpp
      nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpc=tpp
      tpp=tpc/2
   ELSE
      tpp=0
   ENDIF
ENDDO
END SUBROUTINE addtree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine updates the binary tree after the root
! value has been used. The root is replaced by the value
! at the bottom of the tree, which is then filtered down
! to its correct position. This ensures that the tree remains
! balanced.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE downtree
IMPLICIT NONE
INTEGER :: tpp,tpc
REAL(KIND=i10) :: rd1,rd2
TYPE(backpointer) :: exch
!
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
! rd1,rd2 = substitution variables
!
! Replace root of tree with its last value
!
IF(ntr.EQ.1)THEN
   ntr=ntr-1
   RETURN
ENDIF
nsts(btg(ntr)%pz,btg(ntr)%px)=1
btg(1)=btg(ntr)
!
! Reduce size of tree by one
!
ntr=ntr-1
!
! Now filter new root down to its correct position
!
tpp=1
tpc=2*tpp
DO WHILE(tpc.lt.ntr)
!
! Check which of the two children is smallest - use the smallest
!
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
   rd2=ttn(btg(tpc+1)%pz,btg(tpc+1)%px)
   IF(rd1.gt.rd2)THEN
      tpc=tpc+1
   ENDIF
!
!  Check whether the child is smaller than the parent; if so, then swap,
!  if not, then we are done
!
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
   rd2=ttn(btg(tpp)%pz,btg(tpp)%px)
   IF(rd1.lt.rd2)THEN
      nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
      nsts(btg(tpc)%pz,btg(tpc)%px)=tpp
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpp=tpc
      tpc=2*tpp
   ELSE
      tpc=ntr+1
   ENDIF
ENDDO
!
! If ntr is an even number, then we still have one more test to do
!
IF(tpc.eq.ntr)THEN
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
   rd2=ttn(btg(tpp)%pz,btg(tpp)%px)
   IF(rd1.lt.rd2)THEN
      nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
      nsts(btg(tpc)%pz,btg(tpc)%px)=tpp
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
   ENDIF
ENDIF
END SUBROUTINE downtree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine updates a value on the binary tree. The FMM
! should only produce updated values that are less than their
! prior values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE updtree(iz,ix)
IMPLICIT NONE
INTEGER :: ix,iz,tpp,tpc
TYPE(backpointer) :: exch
!
! ix,iz = grid position of new addition to tree
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
!
! Filter the updated value to its correct position
!
tpc=nsts(iz,ix)
tpp=tpc/2
DO WHILE(tpp.gt.0)
   IF(ttn(iz,ix).lt.ttn(btg(tpp)%pz,btg(tpp)%px))THEN
      nsts(iz,ix)=tpp
      nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpc=tpp
      tpp=tpc/2
   ELSE
      tpp=0
   ENDIF
ENDDO
END SUBROUTINE updtree

END MODULE traveltime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to implement the Fast Marching
! Method (FMM) for calculating first-arrival traveltimes
! through a 2-D continuous velocity medium in spherical shell
! coordinates (x=theta or latitude, z=phi or longitude).
! It is written in Fortran 90, although it is probably more
! accurately  described as Fortran 77 with some of the Fortran 90
! extensions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PROGRAM tomo_surf
subroutine CalSurfGAniso(SenGfile,lay,nx,ny,nz,nparpi,vels,iw,rw,col,dsurf, &
               GGc,GGs,GVs,dall,&
              goxdf,gozdf,dvxdf,dvzdf,kmaxRc,tRc,periods,depz,minthk, &
              scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf,nar,writepath,&
              iter,dvsExtral,tRcV)
USE globalp
USE traveltime
IMPLICIT NONE
!CHARACTER (LEN=30) ::grid,frechet
!CHARACTER (LEN=40) :: sources,receivers,otimes
!CHARACTER (LEN=30) :: travelt,rtravel,wrays,cdum
INTEGER :: i,j,k,l,nsrc,tnr,urg
INTEGER :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
INTEGER :: ogx,ogz,grdfx,grdfz,maxbt
REAL(KIND=i10) :: x,z,goxb,gozb,dnxb,dnzb
!--------------------!
! INPUT LIST
! nrcf: maximum number of data for a certain receiver. We output this number when we reformat the data file. You can also set it to nrc*numf, which is the number of stations times the total
!             number of periods for dispersion measurements.
! nsrcsurf=nsrc
! kmax=kmaxRc: total number of period
! nsrcsurf1: <-nsrc1: nsrc1(knum) counter of sources at specific period  (knum: period counter)
! nrc1: nrc1(istep,knum),  counter of records at specific period and source
! rcxf,rczf: rcxf(receiver, source, period) the colatitude, longitude of receiver at fixed source and period ( in rad)
! scxf,sczf: (source, period) the colatitude and longitude of receiver at fixed period ( in rad)
! minthk: 'minthk' = layerthick / sublayernum, can be 1, 2 ,3
! depz: nz, depth grid, load from MOD
! periods: periods(istep,knum) period index for specific source at specific period (istep: at kunm period index, source counter)
! dvxdf,dvzdf: grid interval in lat and lon direction, to build propagation grid
! goxdf,gozdf: goxd gozd (upper left point,[lat,lon])  for inversion grid/ velocity grid
! vels: <-vsf, initial Vs model, (1:nx,1:ny,1:nz) contain edge points.
! nparpi: <-maxvp = (nx-2)*(ny-2)*(nz-1); number of parameter (effective, does not contain edge grid)
! nx ny nz: grid number in lat lon and depth direction
! lay: lay=nz-1
!--------------------!
! OUTPUT LIST
! dsurf: forward predicated data (traveltime-s)
! iw,rw,col: for output G matrix record
! nar
! tRcV: the corresponding phase velocity.
!--------------------!
! sources = File containing source locations
! receivers = File containing receiver locations
! grid = File containing grid of velocity vertices for
!              resampling on a finer grid with cubic B-splines
! frechet = output file containing matrix of frechet derivatives
! travelt = File name for storage of traveltime field
! wttf = Write traveltimes to file? (0=no,>0=source id)
! fom = Use first-order(0) or mixed-order(1) scheme
! nsrc = number of sources
! scx,scz = source location in r,x,z
! scx,scz = source location in r,x,z
! x,z = temporary variables for source location
! fsrt = find source-receiver traveltimes? (0=no,1=yes)
! rtravel = output file for source-receiver traveltimes
! cdum = dummy character variable ! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id.)
! wrays = file containing raypath geometries
! cfd = calculate Frechet derivatives? (0=no, 1=yes)
! tnr = total number of receivers
! sgs = Extent of refined source grid
! isx,isz = cell containing source
! nnxb,nnzb = Backup for nnz,nnx
! goxb,gozb = Backup for gox,goz
! dnxb,dnzb = Backup for dnx,dnz
! ogx,ogz = Location of refined grid origin
! gridfx,grdfz = Number of refined nodes per cell
! urg = use refined grid (0=no,1=yes,2=previously used)
! maxbt = maximum size of narrow band binary tree
! otimes = file containing source-receiver association information
! pvRc= Rayleigh wave phase velocity model, pv(2D point, period)
! velf= 2D phase map at specific period, extracted from pvRc, for the inversion grid
! asgr = Apply source grid refinement (0=no,1=yes), set asgr =1
!c-----------------------------------------------------------------
!	variables defined by Hongjian Fang
	integer nx,ny,nz
       integer kmax,nsrcsurf,nrcf
       real vels(nx,ny,nz)
       real rw(*)
       integer iw(*),col(*)
       real dsurf(*)
       real goxdf,gozdf,dvxdf,dvzdf
       integer kmaxRc
       real*8 tRc(*)
       integer wavetype(nsrcsurf,kmax)
       integer periods(nsrcsurf,kmax),nrc1(nsrcsurf,kmax),nsrcsurf1(kmax)
       integer igrt(nsrcsurf,kmax)
       real scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax),rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
       integer nar
       real minthk
       integer nparpi

	real vpz(nz),vsz(nz),rhoz(nz),depz(nz)
       real*8 pvRc(nx*ny,kmaxRc)
!       real*8 sen_vs(nx*ny,kmax,nz),sen_vp(nx*ny,kmax,nz)
!       real*8 sen_rho(nx*ny,kmax,nz)
!       real coe_rho(nz-1),coe_a(nz-1)
       real*8 velf(ny*nx)
	integer kmax1,kmax2,kmax3,count1
       integer igr
       integer iwave
       integer knumi,srcnum
	real,dimension(:,:),allocatable:: fdm

       real row(3*nparpi)
       real vpft(nz-1)
	real cbst1
        integer ii,jj,kk,nn,istep
       integer level,maxlevel,maxleveld,HorizonType,VerticalType,PorS
      real,parameter::ftol=1e-4
      integer writepath
!    variables defined by Chuanming Liu
        integer dall
       character(len=*)SenGfile
       integer lay
       real*8 Lsen_Vs(nx*ny,kmax,nz-1),Lsen_Gsc(nx*ny,kmax,nz-1)
       real*8 sen_Vs(nx*ny,kmax,nz-1),sen_Gsc(nx*ny,kmax,nz-1)
       real,dimension(:,:),allocatable::fdmc
       real,dimension(:,:),allocatable::fdms
       real rowVs(nparpi),rowGc(nparpi),rowGs(nparpi)
       integer zz,tt
       real GVs(dall,nparpi),GGc(dall,nparpi),GGs(dall,nparpi)

       integer Nperiod
       real*8 Tperiod1,Tperiod2
       logical ex
       integer istat
       character(len=30)Tchar,rayfile

       integer iter
       real dvsExtral(nx-2,ny-2,nz-1)
       real*8 tRcV((nx-2)*(ny-2),kmaxRc)
!--------------------------------------------------------------------------------!
! Part 1. Initialization for assignment of module globap
!--------------------------------------------------------------------------------!
gdx=5
gdz=5
asgr=1
sgdl=8
sgs=8
earth=6371.0
fom=1
snb=0.5
goxd=goxdf
gozd=gozdf
dvxd=dvxdf
dvzd=dvzdf
nvx=nx-2
nvz=ny-2
ALLOCATE(velv(0:nvz+1,0:nvx+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL velv'
ENDIF

!
! Convert from degrees to radians
!
dvx=dvxd*pi/180.0
dvz=dvzd*pi/180.0
gox=(90.0-goxd)*pi/180.0
goz=gozd*pi/180.0
!
! Compute corresponding values for propagation grid.
!
nnx=(nvx-1)*gdx+1
nnz=(nvz-1)*gdz+1
dnx=dvx/gdx
dnz=dvz/gdz
dnxd=dvxd/gdx
dnzd=dvzd/gdz
ALLOCATE(veln(nnz,nnx), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL veln'
ENDIF

!
! Call a subroutine which reads in the velocity grid
!
!CALL gridder(grid)
!
! Read in all source coordinates.
!
!
! Now work out, source by source, the first-arrival traveltime
! field plus source-receiver traveltimes
! and ray paths if required. First, allocate memory to the
! traveltime field array
!
ALLOCATE(ttn(nnz,nnx), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL ttn'
ENDIF
   rbint=0
!
! Allocate memory for node status and binary trees
!
ALLOCATE(nsts(nnz,nnx))
maxbt=NINT(snb*nnx*nnz)
ALLOCATE(btg(maxbt))

!--------------------------------------------------------------------------------!
!  Part 2. Calculation of surface wave sensitivity kernel
!--------------------------------------------------------------------------------!
! For use of layer model to calculate sensitivity kernel,
! nz inversion grid -> nz-1 inversion layer
!  for  old kernel  sen_vsRc(nx*ny,kmaxRc,nz): contain nz grid kernel
! which stands for half space

allocate(fdm(0:nvz+1,0:nvx+1))
allocate(fdmc(0:nvz+1,0:nvx+1))
allocate(fdms(0:nvz+1,0:nvx+1))


open(37,file=SenGfile,action='read',status='old')
DO j=1,nx*ny
   read(37,*)((Lsen_Vs(j,i,zz),zz=1,nz-1),i=1,kmaxRc)
ENDDO
DO j=1,nx*ny
   read(37,*)((Lsen_Gsc(j,i,zz),zz=1,nz-1),i=1,kmaxRc)
ENDDO
close(unit=37)


! cal  phase velocity 3D model.
 iwave=2
 igr=0
 IF (iter==1)THEN
       call CalRayleighPhase(nx,ny,nz,vels,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
 ELSE
        call CalRayPhIteration(nx,ny,nz,vels,dvsExtral,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
ENDIF
!--------------------------------------------------------------------------------!
!  Part 3. Main loop:1.periods-knumi 2.source-srcnum 3.record-istep
!--------------------------------------------------------------------------------!
nar=0
count1=0
sen_Vs=0
sen_Gsc=0
Tperiod1=0

do knumi=1,kmax
do srcnum=1,nsrcsurf1(knumi)
        velf(1:nx*ny)=pvRc(1:nx*ny,periods(srcnum,knumi))
        sen_Vs(:,1:kmaxRc,:)=Lsen_Vs(:,1:kmaxRc,:)!(:,nt(istep),:)
        sen_Gsc(:,1:kmaxRc,:)=Lsen_Gsc(:,1:kmaxRc,:)!(:,nt(istep),:)

call gridder(velf)
   x=scxf(srcnum,knumi)
   z=sczf(srcnum,knumi)
!
!  Begin by computing refined source grid if required
!
   urg=0
   IF(asgr.EQ.1)THEN
!
!     Back up coarse velocity grid to a holding matrix
!
     ALLOCATE(velnb(nnz,nnx))
	! MODIFIEDY BY HONGJIAN FANG @ USTC 2014/04/17
      velnb(1:nnz,1:nnx)=veln(1:nnz,1:nnx)
      nnxb=nnx
      nnzb=nnz
      dnxb=dnx
      dnzb=dnz
      goxb=gox
      gozb=goz
!
!     Identify nearest neighbouring node to source
!
      isx=INT((x-gox)/dnx)+1
      isz=INT((z-goz)/dnz)+1
      sw=0
      IF(isx.lt.1.or.isx.gt.nnx)sw=1
      IF(isz.lt.1.or.isz.gt.nnz)sw=1
      IF(sw.eq.1)then
         x=90.0-x*180.0/pi
         z=z*180.0/pi
         WRITE(6,*)"Source lies outside bounds of model (lat,long)= ",x,z
         WRITE(6,*)"TERMINATING PROGRAM!!!"
         STOP
      ENDIF
      IF(isx.eq.nnx)isx=isx-1
      IF(isz.eq.nnz)isz=isz-1
!
!     Now find rectangular box that extends outward from the nearest source node
!     to "sgs" nodes away.
!
      vnl=isx-sgs
      IF(vnl.lt.1)vnl=1
      vnr=isx+sgs
      IF(vnr.gt.nnx)vnr=nnx
      vnt=isz-sgs
      IF(vnt.lt.1)vnt=1
      vnb=isz+sgs
      IF(vnb.gt.nnz)vnb=nnz
      nrnx=(vnr-vnl)*sgdl+1
      nrnz=(vnb-vnt)*sgdl+1
      drnx=dvx/REAL(gdx*sgdl)
      drnz=dvz/REAL(gdz*sgdl)
      gorx=gox+dnx*(vnl-1)
      gorz=goz+dnz*(vnt-1)
      nnx=nrnx
      nnz=nrnz
      dnx=drnx
      dnz=drnz
      gox=gorx
      goz=gorz
!
!     Reallocate velocity and traveltime arrays if nnx>nnxb or
!     nnz<nnzb.
!
      IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
         idm1=nnx
         IF(nnxb.GT.idm1)idm1=nnxb
         idm2=nnz
         IF(nnzb.GT.idm2)idm2=nnzb
         DEALLOCATE(veln,ttn,nsts,btg)
         ALLOCATE(veln(idm2,idm1))
         ALLOCATE(ttn(idm2,idm1))
         ALLOCATE(nsts(idm2,idm1))
         maxbt=NINT(snb*idm1*idm2)
         ALLOCATE(btg(maxbt))
      ENDIF
!
!     Call a subroutine to compute values of refined velocity nodes
!
      CALL bsplrefine
!
!     Compute first-arrival traveltime field through refined grid.
!
      urg=1
      CALL travel(x,z,urg)
!
!     Now map refined grid onto coarse grid.
!
      ALLOCATE(ttnr(nnzb,nnxb))
      ALLOCATE(nstsr(nnzb,nnxb))
      IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
         idm1=nnx
         IF(nnxb.GT.idm1)idm1=nnxb
         idm2=nnz
         IF(nnzb.GT.idm2)idm2=nnzb
         DEALLOCATE(ttnr,nstsr)
         ALLOCATE(ttnr(idm2,idm1))
         ALLOCATE(nstsr(idm2,idm1))
      ENDIF
      ttnr=ttn
      nstsr=nsts
      ogx=vnl
      ogz=vnt
      grdfx=sgdl
      grdfz=sgdl
      nsts=-1
      DO k=1,nnz,grdfz
         idm1=ogz+(k-1)/grdfz
         DO l=1,nnx,grdfx
            idm2=ogx+(l-1)/grdfx
            nsts(idm1,idm2)=nstsr(k,l)
            IF(nsts(idm1,idm2).GE.0)THEN
               ttn(idm1,idm2)=ttnr(k,l)
            ENDIF
         ENDDO
      ENDDO
!
!     Backup refined grid information
!
      nnxr=nnx
      nnzr=nnz
      goxr=gox
      gozr=goz
      dnxr=dnx
      dnzr=dnz
!
!     Restore remaining values.
!
      nnx=nnxb
      nnz=nnzb
      dnx=dnxb
      dnz=dnzb
      gox=goxb
      goz=gozb
      DO j=1,nnx
         DO k=1,nnz
            veln(k,j)=velnb(k,j)
         ENDDO
      ENDDO
!
!     Ensure that the narrow band is complete; if
!     not, then some alive points will need to be
!     made close.
!
      DO k=1,nnx
         DO l=1,nnz
            IF(nsts(l,k).EQ.0)THEN
               IF(l-1.GE.1)THEN
                  IF(nsts(l-1,k).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(l+1.LE.nnz)THEN
                  IF(nsts(l+1,k).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(k-1.GE.1)THEN
                  IF(nsts(l,k-1).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(k+1.LE.nnx)THEN
                  IF(nsts(l,k+1).EQ.-1)nsts(l,k)=1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Finally, call routine for computing traveltimes once
!     again.
!
      urg=2
      CALL travel(x,z,urg)
   ELSE
!
!     Call a subroutine that works out the first-arrival traveltime
!     field.--- for not Apply source grid refinement
!
      CALL travel(x,z,urg)
   ENDIF
!
!  Find source-receiver traveltimes if required
!
!
   do istep=1,nrc1(srcnum,knumi)
        CALL srtimes(x,z,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),cbst1)
        count1=count1+1
        dsurf(count1)=cbst1
!!-------------------------------------------------------------
!  Calculate raypath geometries and write to file if required.
!  Calculate Frechet derivatives with the same subroutine
!  if required.
! nvz, nvx: inversion grid nvx=nx-2; nvz=ny-2
! fdm(0:nvz+1,0:nvx+1): inversion grid,  raypath sensitivity
! nparpi: (nx-2)*(ny-2)*(nz-1)
! count1: data, d, counter index
! nn: 1: nparpi , model index
! nar: number of G which value is not zero.if Gn*m does not contain zeros,  nar =n*m
! rw: G none zero value, if Gn*m does not contain zeros,  nar =n*mTperiod1
! col: counter index on model parameter, col(i) means model parameter order at i th data in G ---m index
! iw: counter index on data.---n index, iw(i) means data order at i th data in G, i=1: nar
! row: value of  count1 row of G
        Nperiod=periods(srcnum,knumi)
        Tperiod2=tRc(Nperiod)
        IF (writepath.EQ.1)THEN
        IF (abs(Tperiod1-Tperiod2).GT.ftol)then
           inquire(unit=40,exist=ex)
           if(ex) close(unit=40,iostat=istat)
           write(Tchar,'(f5.1)')Tperiod2
           rayfile='raypath_refmdl_'//TRIM(adjustl(Tchar))//'s.dat'
           open(40,file=rayfile,action='write')
           Tperiod1=Tperiod2
         ELSE
            inquire(file=rayfile,exist=ex)
            if(.not.ex) STOP"raypath file hasn't built."
        ENDIF
        ENDIF

        CALL rpathsAzim(x,z,fdm,fdmc,fdms,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),writepath,Tperiod2)
        row(1:3*nparpi)=0.0
        do jj=1,nvz
            do kk=1,nvx
                 if(abs(fdm(jj,kk)).ge.ftol) then
! row for dvs
                   row( (jj-1)*nvx+kk: (nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=sen_Vs(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdm(jj,kk)
! row for Gc/L

                   row(nparpi+ (jj-1)*nvx+kk: nparpi+(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdmc(jj,kk)
! row for Gs/L

                   row(nparpi*2+ (jj-1)*nvx+kk: nparpi*2+(nz-2)*nvz*nvx+(jj-1)*nvx+kk: nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdms(jj,kk)
                endif
            enddo
       enddo

       do nn=1,nparpi*3
	   if(abs(row(nn)).gt.ftol) then
	      nar=nar+1
	      rw(nar)=real(row(nn))
	      iw(nar+1)= count1
              col(nar)=nn
           endif
	enddo
!------------------------------------------------------------------------------------------------------------!
! output G matrix which is used in forward traveltime data.
! NOTE: sen_Gsc:(nx*ny,kmaxRc,nz-1)--> GGcs(dataIndex,nvx*nvz*(nz-1)):
!                                                                               GGcs(dataIndex,2Dgrid-depth1---> depth nz)
        do jj=1,nvz
            do kk=1,nvx
                 if(abs(fdm(jj,kk)).ge.ftol) then
                   GVs(count1, (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Vs(jj*(nvx+2)+kk+1, knumi, 1:nz-1)*fdm(jj,kk)
                   GGc(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdmc(jj,kk)
                   GGs(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdms(jj,kk)
                endif
            enddo
        enddo
!------------------------------------------------------------------------------------------------------------!

   enddo ! records at fixed  source and period.

           IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)

           IF(rbint.EQ.1)THEN
              WRITE(6,*)'Note that at least one two-point ray path'
              WRITE(6,*)'tracked along the boundary of the model.'
              WRITE(6,*)'This class of path is unlikely to be'
              WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
              WRITE(6,*)'that you adjust the dimensions of your grid'
              WRITE(6,*)'to prevent this from occurring.'
           ENDIF
        IF(asgr.EQ.1)THEN
           DEALLOCATE (velnb, STAT=checkstat)
           IF(checkstat > 0)THEN
              WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
           ENDIF
        ENDIF

enddo
enddo

inquire(unit=40,exist=ex)
if(ex) close(unit=40,iostat=istat)
!---------------------------------------------------------------------------------------------------
! output period tomo map in the inversion range
!  pvRc(1:nx*ny,1:kmaxRc)--only need (nvx*nvz))
tRcV=0.0
    DO tt=1,kmaxRc
        DO jj=1,ny-2
            DO ii=1,nx-2
!                WRITE(41,'(5f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),pvRc(jj*nx+ii+1,tt)
                tRcV((jj-1)*(nx-2)+ii,tt)=pvRc(jj*nx+ii+1,tt)
            ENDDO
        ENDDO
    ENDDO

deallocate(fdm)
deallocate(fdmc)
deallocate(fdms)
deallocate(velv,veln,ttn,nsts,btg)
END subroutine



SUBROUTINE gridder(pv)
!subroutine gridder(pv)
!subroutine gridder()
USE globalp
IMPLICIT NONE
INTEGER :: i,j,l,m,i1,j1,conx,conz,stx,stz
REAL(KIND=i10) :: u,sumi,sumj
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi
!CHARACTER (LEN=30) :: grid
!
! u = independent parameter for b-spline
! ui,vi = bspline basis functions
! conx,conz = variables for edge of B-spline grid
! stx,stz = counters for veln grid points
! sumi,sumj = summation variables for computing b-spline
!
!C---------------------------------------------------------------
double precision pv(*)
!integer count1
!C---------------------------------------------------------------
! Open the grid file and read in the velocity grid.
!
!OPEN(UNIT=10,FILE=grid,STATUS='old')
!READ(10,*)nvx,nvz
!READ(10,*)goxd,gozd
!READ(10,*)dvxd,dvzd
!count1=0
DO i=0,nvz+1
   DO j=0,nvx+1
!	count1=count1+1
!      READ(10,*)velv(i,j)
!	velv(i,j)=real(pv(count1))
	velv(i,j)=real(pv(i*(nvx+2)+j+1))
   ENDDO
ENDDO
!CLOSE(10)
!
! Convert from degrees to radians
!
!
! Now dice up the grid
!
ALLOCATE(ui(gdx+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: Subroutine gridder: REAL ui'
ENDIF
DO i=1,gdx+1
   u=gdx
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(gdz+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: Subroutine gridder: REAL vi'
ENDIF
DO i=1,gdz+1
   u=gdz
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
DO i=1,nvz-1
   conz=gdz
   IF(i==nvz-1)conz=gdz+1
   DO j=1,nvx-1
      conx=gdx
      IF(j==nvx-1)conx=gdx+1
      DO l=1,conz
         stz=gdz*(i-1)+l
         DO m=1,conx
            stx=gdx*(j-1)+m
            sumi=0.0
            DO i1=1,4
               sumj=0.0
               DO j1=1,4
                  sumj=sumj+ui(m,j1)*velv(i-2+i1,j-2+j1)
               ENDDO
               sumi=sumi+vi(l,i1)*sumj
            ENDDO
            veln(stz,stx)=sumi
         ENDDO
      ENDDO
   ENDDO
ENDDO
DEALLOCATE(ui,vi, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE gridder: REAL ui,vi'
ENDIF
END SUBROUTINE gridder


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is similar to bsplreg except that it has been
! modified to deal with source grid refinement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bsplrefine
USE globalp
INTEGER :: i,j,k,l,i1,j1,st1,st2,nrzr,nrxr
INTEGER :: origx,origz,conx,conz,idm1,idm2
REAL(KIND=i10) :: u,v
REAL(KIND=i10), DIMENSION (4) :: sum
REAL(KIND=i10), DIMENSION(gdx*sgdl+1,gdz*sgdl+1,4) :: ui,vi
!
! nrxr,nrzr = grid refinement level for source grid in x,z
! origx,origz = local origin of refined source grid
!
! Begin by calculating the values of the basis functions
!
nrxr=gdx*sgdl
nrzr=gdz*sgdl
DO i=1,nrzr+1
   v=nrzr
   v=(i-1)/v
   DO j=1,nrxr+1
      u=nrxr
      u=(j-1)/u
      ui(j,i,1)=(1.0-u)**3/6.0
      ui(j,i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
      ui(j,i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
      ui(j,i,4)=u**3/6.0
      vi(j,i,1)=(1.0-v)**3/6.0
      vi(j,i,2)=(4.0-6.0*v**2+3.0*v**3)/6.0
      vi(j,i,3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
      vi(j,i,4)=v**3/6.0
   ENDDO
ENDDO
!
! Calculate the velocity values.
!
origx=(vnl-1)*sgdl+1
origz=(vnt-1)*sgdl+1
DO i=1,nvz-1
   conz=nrzr
   IF(i==nvz-1)conz=nrzr+1
   DO j=1,nvx-1
      conx=nrxr
      IF(j==nvx-1)conx=nrxr+1
      DO k=1,conz
         st1=gdz*(i-1)+(k-1)/sgdl+1
         IF(st1.LT.vnt.OR.st1.GT.vnb)CYCLE
         st1=nrzr*(i-1)+k
         DO l=1,conx
            st2=gdx*(j-1)+(l-1)/sgdl+1
            IF(st2.LT.vnl.OR.st2.GT.vnr)CYCLE
            st2=nrxr*(j-1)+l
            DO i1=1,4
               sum(i1)=0.0
               DO j1=1,4
                  sum(i1)=sum(i1)+ui(l,k,j1)*velv(i-2+i1,j-2+j1)
               ENDDO
               sum(i1)=vi(l,k,i1)*sum(i1)
            ENDDO
            idm1=st1-origz+1
            idm2=st2-origx+1
            IF(idm1.LT.1.OR.idm1.GT.nnz)CYCLE
            IF(idm2.LT.1.OR.idm2.GT.nnx)CYCLE
            veln(idm1,idm2)=sum(1)+sum(2)+sum(3)+sum(4)
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE bsplrefine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates all receiver traveltimes for
! a given source and writes the results to file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE srtimes(scx,scz,rcx1,rcz1,cbst1)
SUBROUTINE srtimes(scx,scz,rcx1,rcz1,cbst1)
USE globalp
IMPLICIT NONE
INTEGER :: i,k,l,irx,irz,sw,isx,isz,csid
INTEGER, PARAMETER :: noray=0,yesray=1
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(6)
REAL(KIND=i5) :: trr
REAL(KIND=i5), PARAMETER :: norayt=0.0
REAL(KIND=i10) :: drx,drz,produ,scx,scz
REAL(KIND=i10) :: rcx1,rcz1,cbst1
REAL(KIND=i10) :: sred,dpl,rd1,vels,velr
REAL(KIND=i10), DIMENSION (2,2) :: vss
!!------------------------------------------------------
!	modified by Hongjian Fang @ USTC
	integer no_p,nsrc
	real dist
!	real cbst(*) !note that the type difference(kind=i5 vs real)
!	integer cbst_stat(*)
!!------------------------------------------------------
!
! irx,irz = Coordinates of cell containing receiver
! trr = traveltime value at receiver
! produ = dummy multiplier
! drx,drz = receiver distance from (i,j,k) grid node
! scx,scz = source coordinates
! isx,isz = source cell location
! sred = Distance from source to receiver
! dpl = Minimum path length in source neighbourhood.
! vels,velr = velocity at source and receiver
! vss = velocity at four grid points about source or receiver.
! csid = current source ID
! noray = switch to indicate no ray present
! norayt = default value given to null ray
! yesray = switch to indicate that ray is present
!
! Determine source-receiver traveltimes one at a time.
!
!0605DO i=1,nrc
!0605   IF(srs(i,csid).EQ.0)THEN
!0605!      WRITE(10,*)noray,norayt
!0605      CYCLE
!0605   ENDIF
!
!  The first step is to locate the receiver in the grid.
!
   irx=INT((rcx1-gox)/dnx)+1
   irz=INT((rcz1-goz)/dnz)+1
   sw=0
   IF(irx.lt.1.or.irx.gt.nnx)sw=1
   IF(irz.lt.1.or.irz.gt.nnz)sw=1
   IF(sw.eq.1)then
      rcx1=90.0-rcx1*180.0/pi
      rcz1=rcz1*180.0/pi
      WRITE(6,*)"Receiver lies outside model (lat,long)= ",rcx1,rcz1
      WRITE(6,*)"TERMINATING PROGRAM!!!!"
      STOP
   ENDIF
   IF(irx.eq.nnx)irx=irx-1
   IF(irz.eq.nnz)irz=irz-1
!
!  Location of receiver successfully found within the grid. Now approximate
!  traveltime at receiver using bilinear interpolation from four
!  surrounding grid points. Note that bilinear interpolation is a poor
!  approximation when traveltime gradient varies significantly across a cell,
!  particularly near the source. Thus, we use an improved approximation in this
!  case. First, locate current source cell.
!
   isx=INT((scx-gox)/dnx)+1
   isz=INT((scz-goz)/dnz)+1
   dpl=dnx*earth
   rd1=dnz*earth*SIN(gox)
   IF(rd1.LT.dpl)dpl=rd1
   rd1=dnz*earth*SIN(gox+(nnx-1)*dnx)
   IF(rd1.LT.dpl)dpl=rd1
   sred=((scx-rcx1)*earth)**2
   sred=sred+((scz-rcz1)*earth*SIN(rcx1))**2
   sred=SQRT(sred)
   IF(sred.LT.dpl)sw=1
   IF(isx.EQ.irx)THEN
      IF(isz.EQ.irz)sw=1
   ENDIF
   IF(sw.EQ.1)THEN
!
!     Compute velocity at source and receiver
!
      DO k=1,2
         DO l=1,2
            vss(k,l)=veln(isz-1+l,isx-1+k)
         ENDDO
      ENDDO
      drx=(scx-gox)-(isx-1)*dnx
      drz=(scz-goz)-(isz-1)*dnz
      CALL bilinear(vss,drx,drz,vels)
      DO k=1,2
         DO l=1,2
            vss(k,l)=veln(irz-1+l,irx-1+k)
         ENDDO
      ENDDO
      drx=(rcx1-gox)-(irx-1)*dnx
      drz=(rcz1-goz)-(irz-1)*dnz
      CALL bilinear(vss,drx,drz,velr)
      trr=2.0*sred/(vels+velr)
   ELSE
      drx=(rcx1-gox)-(irx-1)*dnx
      drz=(rcz1-goz)-(irz-1)*dnz
      trr=0.0
      DO k=1,2
         DO l=1,2
            produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
            trr=trr+ttn(irz-1+l,irx-1+k)*produ
         ENDDO
      ENDDO
   ENDIF
!   WRITE(10,*)yesray,trr
!!-----------------------------------------------------------------
!	modified bu Hongjian Fang @ USTC
!	count2=count2+1
!	cbst((no_p-1)*nsrc*nrc+(csid-1)*nrc+i)=trr
	cbst1=trr
!	call delsph(scx,scz,rcx(i),rcz(i),dist)
!	travel_path(count2)=dist
!cbst_stat((no_p-1)*nsrc*nrc+(csid-1)*nrc+i)=yesray
!0605ENDDO
END SUBROUTINE srtimes



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed four node values which lie on
! the corners of a rectangle and the coordinates of a point
! lying within the rectangle. It calculates the value at
! the internal point by using bilinear interpolation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bilinear(nv,dsx,dsz,biv)
USE globalp
IMPLICIT NONE
INTEGER :: i,j
REAL(KIND=i10) :: dsx,dsz,biv
REAL(KIND=i10), DIMENSION(2,2) :: nv
REAL(KIND=i10) :: produ
!
! nv = four node vertex values
! dsx,dsz = distance between internal point and top left node
! dnx,dnz = width and height of node rectangle
! biv = value at internal point calculated by bilinear interpolation
! produ = product variable
!
biv=0.0
DO i=1,2
   DO j=1,2
      produ=(1.0-ABS(((i-1)*dnx-dsx)/dnx))*(1.0-ABS(((j-1)*dnz-dsz)/dnz))
      biv=biv+nv(i,j)*produ
   ENDDO
ENDDO
END SUBROUTINE bilinear




subroutine caldespersion(nx,ny,nz,vel,pvRc, &
                  iwave,igr,kmaxRc,tRc,depz,minthk)
        use omp_lib
        implicit none

        integer nx,ny,nz
        real vel(nx,ny,nz)

        integer iwave,igr
        real minthk
        real depz(nz)
        integer kmaxRc
        real*8 tRc(kmaxRc)
        real*8 pvRc(nx*ny,kmaxRc)



        real vpz(nz),vsz(nz),rhoz(nz)
	integer mmax,iflsph,mode,rmax
        integer ii,jj,k,i,nn,kk
    	integer,parameter::NL=200
    	integer,parameter::NP=60
        real*8 cg1(NP),cg2(NP),cga,cgRc(NP)
    	real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    	real depm(NL),vpm(NL),vsm(NL),rhom(NL),thkm(NL)
	real dlnVs,dlnVp,dlnrho


        mmax=nz
        iflsph=1
        mode=1
        dlnVs=0.01
        dlnVp=0.01
        dlnrho=0.01

!$omp parallel &
!$omp default(private) &
!$omp shared(depz,nx,ny,nz,minthk,kmaxRc,mmax,vel) &
!$omp shared(tRc,pvRc,iflsph,iwave,mode,igr)
!$omp do
        do jj=1,ny
    	do ii=1,nx
    	vsz(1:nz)=vel(ii,jj,1:nz)
        ! some other emperical relationship maybe better,
    	do k=1,nz
        vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
        0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
    	rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
               0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
               0.000106*vpz(k)**5
        enddo

    	call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax,rdep,&
        rvp,rvs,rrho,rthk)
    	call surfdisp96(rthk,rvp,rvs,rrho,rmax,iflsph,iwave,mode,igr,kmaxRc,&
        tRc,cgRc)
        pvRc((jj-1)*nx+ii,1:kmaxRc)=cgRc(1:kmaxRc)
        enddo
        enddo
!$omp end do
!$omp end parallel
             end subroutine




! CALSURFGCSG: 1. calculate G dependent matrix for  [dvs/vs0, Gc/L, Gs/L] part correspondingly---
!                                 2. calculate forward traveltime in isotropic true model
! FWDINVRESDATA:
! only target: forward observed traveltime. !!!
! G matrix here not used by next step!
subroutine FwdInvResData(SenGTruefile,lay,nx,ny,nz,nparpi,vels,&
             vsRela,GcInv,GsInv,obsTvs,obsTaa,tRcV,&
              goxdf,gozdf,dvxdf,dvzdf,kmaxRc,tRc,periods,depz,minthk, &
              scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf, dall,writepath)
USE globalp
USE traveltime
IMPLICIT NONE
!CHARACTER (LEN=30) ::grid,frechet
!CHARACTER (LEN=40) :: sources,receivers,otimes
!CHARACTER (LEN=30) :: travelt,rtravel,wrays,cdum
INTEGER :: i,j,k,l,nsrc,tnr,urg
INTEGER :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
INTEGER :: ogx,ogz,grdfx,grdfz,maxbt
REAL(KIND=i10) :: x,z,goxb,gozb,dnxb,dnzb
!--------------------!
! INPUT LIST
! nrcf: maximum number of data for a certain receiver. We output this number when we reformat the data file. You can also set it to nrc*numf, which is the number of stations times the total
!             number of periods for dispersion measurements.
! nsrcsurf=nsrc
! kmax=kmaxRc: total number of period
! nsrcsurf1: <-nsrc1: nsrc1(knum) counter of sources at specific period  (knum: period counter)
! nrc1: nrc1(istep,knum),  counter of records at specific period and source
! rcxf,rczf: rcxf(receiver, source, period) the colatitude, longitude of receiver at fixed source and period ( in rad)
! scxf,sczf: (source, period) the colatitude and longitude of receiver at fixed period ( in rad)
! minthk: 'minthk' = layerthick / sublayernum, can be 1, 2 ,3
! depz: nz, depth grid, load from MOD
! periods: periods(istep,knum) period index for specific source at specific period (istep: at kunm period index, source counter)
! dvxdf,dvzdf: grid interval in lat and lon direction, to build propagation grid
! goxdf,gozdf: goxd gozd (upper left point,[lat,lon])  for inversion grid/ velocity grid
! vels: <-vsf, initial Vs model, (1:nx,1:ny,1:nz) contain edge points.
! nparpi: <-maxvp = (nx-2)*(ny-2)*(nz-1); number of parameter (effective, does not contain edge grid)
! nx ny nz: grid number in lat lon and depth direction
! lay: lay=nz-1
! Add Nov,26
!  Gctrue,Gstrue: true model of Gc, Gs--replaced by GcInv,GsInv
!  dall:  number of data.
! writepath: integer switch of the output raypath file.
!--------------------!
! OUTPUT LIST
! obsTvs: forward predicated traveltime for isotropic model part
! obsTaa: forward predicated traveltime for anisotropy part
! tRcV: forward predicated phase  velocity.
!--------------------!
! sources = File containing source locations
! receivers = File containing receiver locations
! grid = File containing grid of velocity vertices for
!              resampling on a finer grid with cubic B-splines
! frechet = output file containing matrix of frechet derivatives
! travelt = File name for storage of traveltime field
! wttf = Write traveltimes to file? (0=no,>0=source id)
! fom = Use first-order(0) or mixed-order(1) scheme
! nsrc = number of sources
! scx,scz = source location in r,x,z
! scx,scz = source location in r,x,z
! x,z = temporary variables for source location
! fsrt = find source-receiver traveltimes? (0=no,1=yes)
! rtravel = output file for source-receiver traveltimes
! cdum = dummy character variable ! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id.)
! wrays = file containing raypath geometries
! cfd = calculate Frechet derivatives? (0=no, 1=yes)
! tnr = total number of receivers
! sgs = Extent of refined source grid
! isx,isz = cell containing source
! nnxb,nnzb = Backup for nnz,nnx
! goxb,gozb = Backup for gox,goz
! dnxb,dnzb = Backup for dnx,dnz
! ogx,ogz = Location of refined grid origin
! gridfx,grdfz = Number of refined nodes per cell
! urg = use refined grid (0=no,1=yes,2=previously used)
! maxbt = maximum size of narrow band binary tree
! otimes = file containing source-receiver association information
! pvRc= Rayleigh wave phase velocity model, pv(2D point, period)
! velf= 2D phase map at specific period, extracted from pvRc, for the inversion grid
! asgr = Apply source grid refinement (0=no,1=yes), set asgr =1
! Added in Nov,24
! ! GVs,GGc,GGs: indirect m*n G matrix for [dvs/vs0, Gc/L, Gs/L]
!c-----------------------------------------------------------------
!	variables defined by Hongjian Fang
	integer nx,ny,nz
       integer kmax,nsrcsurf,nrcf
       real vels(nx,ny,nz)
!       real rw(*)
!       integer iw(*),col(*)
       real goxdf,gozdf,dvxdf,dvzdf
       integer kmaxRc
       real*8 tRc(*)
       integer wavetype(nsrcsurf,kmax)
       integer periods(nsrcsurf,kmax),nrc1(nsrcsurf,kmax),nsrcsurf1(kmax)
       integer igrt(nsrcsurf,kmax)
       real scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax),rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
!       integer nar
       real minthk
       integer nparpi

	real vpz(nz),vsz(nz),rhoz(nz),depz(nz)
       real*8 pvRc(nx*ny,kmaxRc)

!       real*8 sen_vs(nx*ny,kmax,nz),sen_vp(nx*ny,kmax,nz)
!       real*8 sen_rho(nx*ny,kmax,nz)
!       real coe_rho(nz-1),coe_a(nz-1)
       real*8 velf(ny*nx)
	integer kmax1,kmax2,kmax3,count1
       integer igr
       integer iwave
       integer knumi,srcnum
	real,dimension(:,:),allocatable:: fdm

       real row(3*nparpi)
       real vpft(nz-1)
	real cbst1
        integer ii,jj,kk,nn,istep
       integer level,maxlevel,maxleveld,HorizonType,VerticalType,PorS
      real,parameter::ftol=1e-4
      integer writepath
!    variables defined by Chuanming Liu
       real GcInv(nx-2,ny-2,nz-1),GsInv(nx-2,ny-2,nz-1)
       integer dall
       character(len=*)SenGTruefile
       integer lay
       real*8 Lsen_Vs(nx*ny,kmax,nz-1),Lsen_Gsc(nx*ny,kmax,nz-1)
       real*8 sen_Vs(nx*ny,kmax,nz-1),sen_Gsc(nx*ny,kmax,nz-1)
       real,dimension(:,:),allocatable::fdmc
       real,dimension(:,:),allocatable::fdms
       real rowVs(nparpi),rowGc(nparpi),rowGs(nparpi)
       integer zz,tt
       real obsTvs(*)
       real obsTaa(*)
       real obsTgc(dall),obsTgs(dall), obst(dall)
       real GVs(dall,nparpi),GGc(dall,nparpi),GGs(dall,nparpi)
       real GcCol(nparpi),GsCol(nparpi)

       integer Nperiod
       real*8 Tperiod1,Tperiod2
       logical ex
       integer istat
       character(len=30)Tchar,rayfile
       real*8 tRcV((nx-2)*(ny-2),kmaxRc)
       real vsRela(nx-2,ny-2,nz-1)
!--------------------------------------------------------------------------------!
! Part 1. Initialization for assignment of module globap
!--------------------------------------------------------------------------------!

gdx=5
gdz=5
asgr=1
sgdl=8
sgs=8
earth=6371.0
fom=1
snb=0.5
goxd=goxdf
gozd=gozdf
dvxd=dvxdf
dvzd=dvzdf
nvx=nx-2
nvz=ny-2
ALLOCATE(velv(0:nvz+1,0:nvx+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL velv'
ENDIF
!
! Convert from degrees to radians
!
dvx=dvxd*pi/180.0
dvz=dvzd*pi/180.0
gox=(90.0-goxd)*pi/180.0
goz=gozd*pi/180.0
!
! Compute corresponding values for propagation grid.
!
nnx=(nvx-1)*gdx+1
nnz=(nvz-1)*gdz+1
dnx=dvx/gdx
dnz=dvz/gdz
dnxd=dvxd/gdx
dnzd=dvzd/gdz
ALLOCATE(veln(nnz,nnx), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL veln'
ENDIF

!
! Call a subroutine which reads in the velocity grid
!
!CALL gridder(grid)
!
! Read in all source coordinates.
!
!
! Now work out, source by source, the first-arrival traveltime
! field plus source-receiver traveltimes
! and ray paths if required. First, allocate memory to the
! traveltime field array
!
ALLOCATE(ttn(nnz,nnx), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL ttn'
ENDIF
   rbint=0
!
! Allocate memory for node status and binary trees
!
ALLOCATE(nsts(nnz,nnx))
maxbt=NINT(snb*nnx*nnz)
ALLOCATE(btg(maxbt))

!--------------------------------------------------------------------------------!
!  Part 2. Calculation of surface wave sensitivity kernel
!--------------------------------------------------------------------------------!
! For use of layer model to calculate sensitivity kernel,
! nz inversion grid -> nz-1 inversion layer
!  for  old kernel  sen_vsRc(nx*ny,kmaxRc,nz): contain nz grid kernel
! which stands for half space

allocate(fdm(0:nvz+1,0:nvx+1))
allocate(fdmc(0:nvz+1,0:nvx+1))
allocate(fdms(0:nvz+1,0:nvx+1))


open(37,file=SenGTruefile,action='read',status='old')
DO j=1,nx*ny
   read(37,*)((Lsen_Vs(j,i,zz),zz=1,nz-1),i=1,kmaxRc)
ENDDO
DO j=1,nx*ny
   read(37,*)((Lsen_Gsc(j,i,zz),zz=1,nz-1),i=1,kmaxRc)
ENDDO
close(unit=37)

! cal  phase velocity 3D model.
 iwave=2
 igr=0
! call CalRayleighPhase(nx,ny,nz,vels,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
 call CalRayPhIteration(nx,ny,nz,vels,vsRela,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
!--------------------------------------------------------------------------------!
!  Part 3. Main loop:1.periods-knumi 2.source-srcnum 3.record-istep
!--------------------------------------------------------------------------------!
count1=0
sen_Vs=0
sen_Gsc=0
Tperiod1=0

do knumi=1,kmax
do srcnum=1,nsrcsurf1(knumi)
        velf(1:nx*ny)=pvRc(1:nx*ny,periods(srcnum,knumi))
        sen_Vs(:,1:kmaxRc,:)=Lsen_Vs(:,1:kmaxRc,:)!(:,nt(istep),:)
        sen_Gsc(:,1:kmaxRc,:)=Lsen_Gsc(:,1:kmaxRc,:)!(:,nt(istep),:)

call gridder(velf)
   x=scxf(srcnum,knumi)
   z=sczf(srcnum,knumi)
!
!  Begin by computing refined source grid if required
!
   urg=0
   IF(asgr.EQ.1)THEN
!
!     Back up coarse velocity grid to a holding matrix
!
     ALLOCATE(velnb(nnz,nnx))
	! MODIFIEDY BY HONGJIAN FANG @ USTC 2014/04/17
      velnb(1:nnz,1:nnx)=veln(1:nnz,1:nnx)
      nnxb=nnx
      nnzb=nnz
      dnxb=dnx
      dnzb=dnz
      goxb=gox
      gozb=goz
!
!     Identify nearest neighbouring node to source
!
      isx=INT((x-gox)/dnx)+1
      isz=INT((z-goz)/dnz)+1
      sw=0
      IF(isx.lt.1.or.isx.gt.nnx)sw=1
      IF(isz.lt.1.or.isz.gt.nnz)sw=1
      IF(sw.eq.1)then
         x=90.0-x*180.0/pi
         z=z*180.0/pi
         WRITE(6,*)"Source lies outside bounds of model (lat,long)= ",x,z
         WRITE(6,*)"TERMINATING PROGRAM!!!"
         STOP
      ENDIF
      IF(isx.eq.nnx)isx=isx-1
      IF(isz.eq.nnz)isz=isz-1
!
!     Now find rectangular box that extends outward from the nearest source node
!     to "sgs" nodes away.
!
      vnl=isx-sgs
      IF(vnl.lt.1)vnl=1
      vnr=isx+sgs
      IF(vnr.gt.nnx)vnr=nnx
      vnt=isz-sgs
      IF(vnt.lt.1)vnt=1
      vnb=isz+sgs
      IF(vnb.gt.nnz)vnb=nnz
      nrnx=(vnr-vnl)*sgdl+1
      nrnz=(vnb-vnt)*sgdl+1
      drnx=dvx/REAL(gdx*sgdl)
      drnz=dvz/REAL(gdz*sgdl)
      gorx=gox+dnx*(vnl-1)
      gorz=goz+dnz*(vnt-1)
      nnx=nrnx
      nnz=nrnz
      dnx=drnx
      dnz=drnz
      gox=gorx
      goz=gorz
!
!     Reallocate velocity and traveltime arrays if nnx>nnxb or
!     nnz<nnzb.
!
      IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
         idm1=nnx
         IF(nnxb.GT.idm1)idm1=nnxb
         idm2=nnz
         IF(nnzb.GT.idm2)idm2=nnzb
         DEALLOCATE(veln,ttn,nsts,btg)
         ALLOCATE(veln(idm2,idm1))
         ALLOCATE(ttn(idm2,idm1))
         ALLOCATE(nsts(idm2,idm1))
         maxbt=NINT(snb*idm1*idm2)
         ALLOCATE(btg(maxbt))
      ENDIF
!
!     Call a subroutine to compute values of refined velocity nodes
!
      CALL bsplrefine
!
!     Compute first-arrival traveltime field through refined grid.
!
      urg=1
      CALL travel(x,z,urg)
!
!     Now map refined grid onto coarse grid.
!
      ALLOCATE(ttnr(nnzb,nnxb))
      ALLOCATE(nstsr(nnzb,nnxb))
      IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
         idm1=nnx
         IF(nnxb.GT.idm1)idm1=nnxb
         idm2=nnz
         IF(nnzb.GT.idm2)idm2=nnzb
         DEALLOCATE(ttnr,nstsr)
         ALLOCATE(ttnr(idm2,idm1))
         ALLOCATE(nstsr(idm2,idm1))
      ENDIF
      ttnr=ttn
      nstsr=nsts
      ogx=vnl
      ogz=vnt
      grdfx=sgdl
      grdfz=sgdl
      nsts=-1
      DO k=1,nnz,grdfz
         idm1=ogz+(k-1)/grdfz
         DO l=1,nnx,grdfx
            idm2=ogx+(l-1)/grdfx
            nsts(idm1,idm2)=nstsr(k,l)
            IF(nsts(idm1,idm2).GE.0)THEN
               ttn(idm1,idm2)=ttnr(k,l)
            ENDIF
         ENDDO
      ENDDO
!
!     Backup refined grid information
!
      nnxr=nnx
      nnzr=nnz
      goxr=gox
      gozr=goz
      dnxr=dnx
      dnzr=dnz
!
!     Restore remaining values.
!
      nnx=nnxb
      nnz=nnzb
      dnx=dnxb
      dnz=dnzb
      gox=goxb
      goz=gozb
      DO j=1,nnx
         DO k=1,nnz
            veln(k,j)=velnb(k,j)
         ENDDO
      ENDDO
!
!     Ensure that the narrow band is complete; if
!     not, then some alive points will need to be
!     made close.
!
      DO k=1,nnx
         DO l=1,nnz
            IF(nsts(l,k).EQ.0)THEN
               IF(l-1.GE.1)THEN
                  IF(nsts(l-1,k).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(l+1.LE.nnz)THEN
                  IF(nsts(l+1,k).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(k-1.GE.1)THEN
                  IF(nsts(l,k-1).EQ.-1)nsts(l,k)=1
               ENDIF
               IF(k+1.LE.nnx)THEN
                  IF(nsts(l,k+1).EQ.-1)nsts(l,k)=1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Finally, call routine for computing traveltimes once
!     again.
!
      urg=2
      CALL travel(x,z,urg)
   ELSE
!
!     Call a subroutine that works out the first-arrival traveltime
!     field.--- for not Apply source grid refinement
!
      CALL travel(x,z,urg)
   ENDIF
!
!  Find source-receiver traveltimes if required
!
!

   do istep=1,nrc1(srcnum,knumi)
        CALL srtimes(x,z,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),cbst1)
        count1=count1+1
        obsTvs(count1)=cbst1
!!-------------------------------------------------------------
!  Calculate raypath geometries and write to file if required.
!  Calculate Frechet derivatives with the same subroutine
!  if required.
! nvz, nvx: inversion grid nvx=nx-2; nvz=ny-2
! fdm(0:nvz+1,0:nvx+1): inversion grid,  raypath sensitivity
! nparpi: (nx-2)*(ny-2)*(nz-1)
! count1: data, d, counter index
! nn: 1: nparpi , model index
! nar: number of G which value is not zero.if Gn*m does not contain zeros,  nar =n*m
! rw: G none zero value, if Gn*m does not contain zeros,  nar =n*m
! col: counter index on model parameter, col(i) means model parameter order at i th data in G ---m index
! iw: counter index on data.---n index, iw(i) means data order at i th data in G, i=1: nar
        Nperiod=periods(srcnum,knumi)
        Tperiod2=tRc(Nperiod)
        IF (writepath.EQ.1)THEN
        IF (abs(Tperiod1-Tperiod2).GT.ftol)then
           inquire(unit=40,exist=ex)
           if(ex) close(unit=40,iostat=istat)
           write(Tchar,'(f5.1)')Tperiod2
           rayfile='raypath_Truemdl_'//TRIM(adjustl(Tchar))//'s.dat'
           open(40,file=rayfile,action='write')
           Tperiod1=Tperiod2
        ELSE
            inquire(file=rayfile,exist=ex)
            if(.not.ex) STOP"raypath file hasn't built."
        ENDIF
        ENDIF
        CALL rpathsAzim(x,z,fdm,fdmc,fdms,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),writepath,Tperiod2)
! NOTE: sen_Gsc:(nx*ny,kmaxRc,nz-1)--> GGcs(dataIndex,nvx*nvz*(nz-2))
        do jj=1,nvz
            do kk=1,nvx
                 if(abs(fdm(jj,kk)).ge.ftol) then
                   GVs(count1, (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Vs(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdm(jj,kk)
                   GGc(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdmc(jj,kk)
                   GGs(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   sen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdms(jj,kk)
                endif
            enddo
        enddo


   enddo ! records at fixed  source and period.

   IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)

   IF(rbint.EQ.1)THEN
       WRITE(6,*)'Note that at least one two-point ray path'
       WRITE(6,*)'tracked along the boundary of the model.'
       WRITE(6,*)'This class of path is unlikely to be'
       WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
       WRITE(6,*)'that you adjust the dimensions of your grid'
       WRITE(6,*)'to prevent this from occurring.'
   ENDIF
   IF(asgr.EQ.1)THEN
       DEALLOCATE (velnb, STAT=checkstat)
       IF(checkstat > 0)THEN
            WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
       ENDIF
   ENDIF
enddo
enddo
!---------------------------------------------------------------------------------------------------

inquire(unit=40,exist=ex)
if(ex) close(unit=40,iostat=istat)
! forward data
! GGc: n*m n: dall   m=maxvp
! GcCol: 1*m:
! Gm=t_ani

!nvx=nx-2
!nvz=ny-2

DO jj=1,ny-2
  DO kk=1,nx-2
       GcCol( (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=GcInv(kk,jj,1:nz-1)
       GsCol( (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=GsInv(kk,jj,1:nz-1)
  ENDDO
ENDDO

obsTgc =MATMUL(GGc,GcCol)
obsTgs=MATMUL(GGs,GsCol)

DO i=1,dall
 obsTaa(i)=obsTgc(i)+obsTgs(i)
 obst(i)=obsTvs(i)+ obsTaa(i)
ENDDO
!---------------------------------------------------------------------------------------------------
! output period tomo map in the inversion range
!  pvRc(1:nx*ny,1:kmaxRc)--only need (nvx*nvz))
tRcV=0.0
    DO tt=1,kmaxRc
        DO jj=1,ny-2
            DO ii=1,nx-2
!                WRITE(41,'(5f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),pvRc(jj*nx+ii+1,tt)
                tRcV((jj-1)*(nx-2)+ii,tt)=pvRc(jj*nx+ii+1,tt)
            ENDDO
        ENDDO
    ENDDO

deallocate(fdm)
deallocate(fdmc)
deallocate(fdms)
deallocate(velv,veln,ttn,nsts,btg)
END subroutine
