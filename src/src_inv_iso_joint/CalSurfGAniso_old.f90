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

   end subroutine


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


  subroutine SpliceRefinedKernel(nx,ny,nz,rmax,vel,kmaxRc,depz,minthk,Lsen_dcdL,Lsen_dcdA,Lsen_Gsc)
        use omp_lib
        implicit none
! INPUT
        integer nx,ny,nz
        real vel(nx,ny,nz)
        integer rmax,rmax1
        integer iwave,igr
        real minthk
        real depz(nz)
        integer kmaxRc
! OUTPUT
        real*4 Lsen_dcdA(nx*ny,kmaxRc,rmax),Lsen_dcdL(nx*ny,kmaxRc,rmax)
        real*4 Lsen_Gsc(nx*ny,kmaxRc,nz-1)
!            real*4,dimension(:,:,:),allocatable::Lsen_dcdL,Lsen_dcdA,Lsen_Gsc


        REAL*4 dcR_dL(kmaxRc,rmax),dcR_dA(kmaxRc,rmax)

! PARAMETETER LIST
        real vpz(nz),vsz(nz),rhoz(nz)
	integer mmax,iflsph,mode
        integer ii,jj,k,i,j,nn,kk,jjj
    	integer,parameter::NL=200
    	real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    	real depm(NL),vpm(NL),vsm(NL),rhom(NL),thkm(NL)
	real dlnVs,dlnVp,dlnrho
        integer nsublay(NL)
         real*4:: A_TIM(rmax),L_TIM(rmax)
        integer post,checkstat

!allocate(Lsen_dcdL(nx*ny,kmaxRc,rmax),Lsen_dcdA(nx*ny,kmaxRc,rmax),stat=checkstat)
!allocate(Lsen_Gsc(nx*ny,kmaxRc,nz-1),stat=checkstat)

        mmax=nz
         Lsen_Gsc=0.0
        write(*,*) 'First Satance Start!'
        write(*,*)'rmax=',rmax
!$omp parallel &
!$omp default(private) &
!$omp shared(depz,nx,ny,nz,minthk,kmaxRc,rmax,mmax,vel) &
!$omp shared(Lsen_dcdL,Lsen_dcdA,Lsen_Gsc)
!$omp do
        DO jj=1,ny
    	   DO ii=1,nx
    	       post=ii+(jj-1)*nx
    	       dcR_dL(1:kmaxRc,1:rmax)=Lsen_dcdL(post,1:kmaxRc,1:rmax)
    	       dcR_dA(1:kmaxRc,1:rmax)=Lsen_dcdA(post,1:kmaxRc,1:rmax)

    	       vsz(1:nz)=vel(ii,jj,1:nz)
        ! some other emperical relationship maybe better,
               DO k=1,nz
                   vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
                   0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
    	           rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
                   0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
                   0.000106*vpz(k)**5
               ENDDO

               CALL refineLayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax1,rdep,&
               rvp,rvs,rrho,rthk, nsublay)
               IF (rmax.NE.rmax1) STOP 'Different rmax: Subroutine CalSenKernel.'
               DO i=1,kmaxRc  ! period
                   k=0
                   DO j=1,nz-1     ! inversion layer
                       DO jjj=1,nsublay(j)    ! refined layer jj th in jth inversion layer
                           k=k+1
                           A_TIM(k)=rrho(k)*rvp(k)*rvp(k)
                           L_TIM(k)=rrho(k)*rvs(k)*rvs(k)
                           Lsen_Gsc(post,i,j)=Lsen_Gsc(post,i,j)+dcR_dA(i,k)*A_TIM(k)+dcR_dL(i,k)*L_TIM(k)
                       ENDDO
                   ENDDO
               ENDDO


           ENDDO
    	ENDDO
!$omp end do
!$omp end parallel

!deallocate(Lsen_dcdL,Lsen_dcdA,Lsen_Gsc)
   end subroutine



subroutine refineLayerMdl(minthk0,mmax,dep,vp,vs,rho,&
                  rmax,rdep,rvp,rvs,rrho,rthk,nsublay)
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
!        k = k + 1
!        rthk(k) = 0.0
!        rvp(k) = vp(mmax)
!        rvs(k) = vs(mmax)
!        rrho(k) = rho(mmax)
!	rdep(k) = dep(mmax)

        rmax = k

        return
        end subroutine


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
subroutine CalSurfGAniso(SenGfile,nx,ny,nz,nparpi,vels,iw,rw,col,dsurf, &
               GGc,GGs,Lsen_Gsc,dall,rmax,tRcV,&
              goxdf,gozdf,dvxdf,dvzdf,kmaxRc,tRc,periods,depz,minthk, &
              scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf,nar,writepath)
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

       real row(2*nparpi)
       real vpft(nz-1)
	real cbst1
        integer ii,jj,kk,nn,istep
       integer level,maxlevel,maxleveld,HorizonType,VerticalType,PorS
      real,parameter::ftol=1e-4
      integer writepath
!    variables defined by Chuanming Liu
        integer dall
       character(len=*)SenGfile
       integer lay,rmax
!       real*8 Lsen_dcdL(nx*ny,kmax,rmax),Lsen_dcdA(nx*ny,kmax,rmax)
!       real*8 sen_Gsc(nx*ny,kmax,nz-1),Lsen_Gsc(nx*ny,kmax,nz-1)

       real*4 Lsen_dcdL(nx*ny,kmax,rmax),Lsen_dcdA(nx*ny,kmax,rmax),Lsen_Gsc(nx*ny,kmax,nz-1)


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
   read(37,*)((Lsen_dcdL(j,i,zz),zz=1,rmax),i=1,kmaxRc)
ENDDO
DO j=1,nx*ny
   read(37,*)((Lsen_dcdA(j,i,zz),zz=1,rmax),i=1,kmaxRc)
ENDDO
close(unit=37)
write(6,*)'Load Sensitivity kernel file successfully!'

Lsen_Gsc=0
write(6,*)'SpliceRefinedKernel Begin!'
call SpliceRefinedKernel(nx,ny,nz,rmax,vels,kmaxRc,depz,minthk,Lsen_dcdL,Lsen_dcdA,Lsen_Gsc)
write(6,*)'SpliceRefinedKernel successfully!'

! cal  phase velocity 3D model.
 iwave=2
 igr=0
 call CalRayleighPhase(nx,ny,nz,vels,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk)
 write(6,*)'CalRayleighPhase successfully!'

! !-----------------------------------------------------!
!            open(37,file='TEST_SenGsc.dat',action='write')
!            DO j=1,nx*ny
!               write(37,*)((Lsen_Gsc(j,i,zz),zz=1,nz-1),i=1,kmaxRc)
!            ENDDO
!            close(37)

!------------------------------------------------------!



!--------------------------------------------------------------------------------!
!  Part 3. Main loop:1.periods-knumi 2.source-srcnum 3.record-istep
!--------------------------------------------------------------------------------!
nar=0
count1=0
!sen_Gsc=0
Tperiod1=0
 write(*,*)'Begin cycle on data!'

do knumi=1,kmax

do srcnum=1,nsrcsurf1(knumi)

        velf(1:nx*ny)=pvRc(1:nx*ny,periods(srcnum,knumi))

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
            if(.not.ex) STOP "raypath file hasn't built."
        ENDIF
        ENDIF

        CALL rpathsAzim(x,z,fdm,fdmc,fdms,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),writepath,Tperiod2)
        row(1:2*nparpi)=0.0
        do jj=1,nvz
            do kk=1,nvx
                 if(abs(fdm(jj,kk)).ge.ftol) then
! row for Gc/L
                   row( (jj-1)*nvx+kk: (nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   Lsen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdmc(jj,kk)
! row for Gs/L
                   row(nparpi+ (jj-1)*nvx+kk: nparpi+(nz-2)*nvz*nvx+(jj-1)*nvx+kk: nvx*nvz)=&
                   Lsen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdms(jj,kk)
                endif
            enddo
       enddo

       do nn=1,nparpi*2
	   if(abs(row(nn)).gt.ftol) then
	      nar=nar+1
	      rw(nar)=real(row(nn))
	      iw(nar+1)= count1
              col(nar)=nn
           else
               row(nn)=0.0
           endif
	enddo
!-----------------------------------------------------!
! FOR TEST
!               write(37,*)(row(i),i=1,2*nparpi)
!------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------!
! output G matrix which is used in forward Calculate ||W(Gm-d)||2
! NOTE: sen_Gsc:(nx*ny,kmaxRc,nz-1)--> GGcs(dataIndex,nvx*nvz*(nz-1)):
!                                                                               GGcs(dataIndex,2Dgrid-depth1---> depth nz)
        do jj=1,nvz
            do kk=1,nvx
                 if(abs(fdm(jj,kk)).ge.ftol) then

                   GGc(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   Lsen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdmc(jj,kk)

                   GGs(count1,(jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                   Lsen_Gsc(jj*(nvx+2)+kk+1,knumi,1:nz-1)*fdms(jj,kk)
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

inquire(unit=37,exist=ex)
if(ex) close(unit=37,iostat=istat)
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
