! Direct Inversion for 3-D Azimuthal Anisotropy
! Forward output part
! V1.1 (2017) Based on Mineos to calculate kernel
! V1.2 (2019) Based on CPS-tregn96 to calcualte kernel
! Copyright:
!    Author: Chuanming Liu (at CU Boulder)
!     Email: Chuanming.liu@colorado.edu
! Reference:
program SurfAnisoForward
    use omp_lib
    implicit none

    character inputfile*80
    character logfile*100
    character outmodel*100
    character outsyn*100
    logical ex
    character dummy*40
    character datafile*80

    integer nx,ny,nz
    real goxd,gozd
    real dvxd,dvzd
    integer nsrc,nrc
    real weight,weight0
    real damp
    real minthk
    integer kmax,kmaxRc
    real*8,dimension(:),allocatable:: tRc
    real,dimension(:),allocatable:: depz
    integer itn
    integer nout
    integer localSize
    real mean,std_devs,balances,balanceb
    integer msurf
    real,dimension(:),allocatable:: obst,dsyn,cbst,dist,periodRre
    real,dimension(:),allocatable:: pvall
    real sta1_lat,sta1_lon,sta2_lat,sta2_lon
    real dist1
    integer dall
    integer istep
    real,parameter :: pi=3.1415926535898
    integer checkstat
    integer ii,jj,kk
    real, dimension (:,:), allocatable :: scxf,sczf
    real, dimension (:,:,:), allocatable :: rcxf,rczf
    integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
    integer,dimension(:),allocatable::nsrc1
    integer,dimension(:,:),allocatable::periods
    real,dimension(:),allocatable::rw
    integer,dimension(:),allocatable::iw,col
    real,dimension(:),allocatable::dv,norm
    real,dimension(:,:,:),allocatable::vsf
    real,dimension(:,:,:),allocatable::vsftrue
    character strf
    integer veltp,wavetp
    real velvalue
    integer knum,knumo,err
    integer istep1,istep2
    integer period
    integer knumi,srcnum,count1
    integer HorizonType,VerticalType
    character line*200
    integer iter,maxiter
    integer maxnar
    real acond
    real anorm
    real arnorm
    real rnorm
    real xnorm
    character str1
    real atol,btol
    real conlim
    integer istop
    integer itnlim
    integer lenrw,leniw
    integer nar,nar_tmp,nars
    integer count3,nvz,nvx
    integer m,maxvp,n
    integer i,j,k
    real spfra
    real noiselevel
    integer ifsyn

    real averdws
    real maxnorm
    !	    real threshold,threshold0
    ! Added by C.M. L
    character SenGTruefile*200, SenGfile*200
    integer lay
    integer maxm
    !            real, dimension (:), allocatable :: deltaT
    real, dimension(:), allocatable:: sigmaT
    real, dimension(:), allocatable:: Tdata,resbst
    real meandeltaT,stddeltaT
    real guassVs,LcorrXY,sigmaDeltaVs,sigmaGcs
    real LcorrZZ
    real,dimension(:,:,:),allocatable:: Gctrue,Gstrue
    real,dimension(:,:),allocatable:: GVs,GGs,GGc
    real,dimension(:,:,:), allocatable :: Lsen_Gsc
    real,dimension(:),allocatable:: obsTvs,obsTaa
    real gaussian
    external gaussian
    real,dimension(:),allocatable:: synT,noise
    real,dimension(:,:,:),allocatable::vsRela,vsAbs,vsReal
    real,dimension(:,:,:),allocatable::gcf,gsf
    integer nar1
    real vsref
    real res2Norm,Mnorm2
    integer,dimension(:),allocatable::iw2
    integer writeperiod
    real*8,dimension(:,:),allocatable:: tRcV
    real sumObs,sumNoise,sumAdd
    logical::  writepath
    real,dimension(:),allocatable:: obst0
    real sta1_latD,sta1_lonD
    real Tvalue,velTrue,sumTnos
    integer rmax
    real*8 :: startT,endT


    startT=OMP_get_wtime()
    write(*,*)
    write(*,*),'                SurfAniso Forward'

    ! load contral parameter file
    if (iargc() < 1) then
        write(*,*) 'input file [SurfAniso.in(default)]:'
        read(*,'(a)') inputfile
        if (len_trim(inputfile) <=1 ) then
            inputfile = 'SurfAnisoForward.in'
        else
            inputfile = inputfile(1:len_trim(inputfile))
        endif
    else
        call getarg(1,inputfile)
    endif

    inquire(file = inputfile, exist = ex)
    if (.not. ex) stop 'unable to open the inputfile'

    open(10,file=inputfile, status='old', action='read')
    read(10,'(a30)')dummy
    read(10,'(a30)')dummy
    read(10,'(a30)')dummy
    read(10,*) datafile
    read(10,*) nx,ny,nz
    read(10,*) goxd,gozd
    read(10,*) dvxd,dvzd
    read(10,*) nsrc
    read(10,*) minthk
    read(10,*) spfra
    read(10,*) writepath
    read(10,*) kmaxRc

    write(*,*)'input Rayleigh wave phase velocity data file:'
    write(*,'(a)') datafile
    write(*,*)  'model origin:latitude,longitue'
    write(*,'(2f10.4)') goxd, gozd
    write(*,*) 'grid spacing:latitude,longitue'
    write(*,'(2f10.4)') dvxd, dvzd
    write(*,*) 'model dimension:nx,ny,nz'
    write(*,'(3i5)') nx, ny, nz
    write(*,*)'depth refined interval layer '
    write(*,'(f8.1)') minthk
    write(*,*) 'number of period'
    write(*,'(i6)') kmaxRc

	write(logfile,'(a,a)')trim(inputfile),'.log'
    open(66, file=logfile,action='write')
    write(66,*)
    write(66,*),'                    SurfAnisoForward'
    write(66,*) 'model origin:latitude,longitue'
    write(66,'(2f10.4)') goxd,gozd
    write(66,*) 'grid spacing:latitude,longitue'
    write(66,'(2f10.4)') dvxd,dvzd
    write(66,*) 'model dimension:nx,ny,nz'
    write(66,'(3i5)') nx,ny,nz

    if (kmaxRc.gt.0) then
        allocate(tRc(kmaxRc), stat=checkstat)
        if (checkstat > 0) stop 'error allocating RP'
        read(10,*)(tRc(i),i=1,kmaxRc)
        write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
        write(*,'(50f6.2)')(tRc(i),i=1,kmaxRc)
        write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
        write(66,'(50f6.2)')(tRc(i),i=1,kmaxRc)
    else
        stop 'Can only deal with Rayleigh wave phase velocity data!'
    endif
    nrc=nsrc
    kmax=kmaxRc
    lay=nz-1
    !--------------------------------------------------------------------------!
    ! READ MEASUREMENTS
    inquire(file=datafile, exist=ex)
    if (.not. ex) stop 'unable to open the datafile'
    write(*,*)'begin load data file.....'

    open(unit=87, file=datafile,status='old')
    allocate(scxf(nsrc,kmax),sczf(nsrc,kmax),&
    rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),stat=checkstat)
    if(checkstat > 0)then
        write(6,*)'scxf error with allocate'
        write(6,*)'nsrc=',nsrc,'kmax=',kmax
    endif

    allocate(periods(nsrc,kmax),wavetype(nsrc,kmax),&
        nrc1(nsrc,kmax),nsrc1(kmax),&
        igrt(nsrc,kmax),stat=checkstat)
    if(checkstat > 0)then
        write(6,*)'error with allocate'
    endif
    allocate(obst(nrc*nsrc*kmax),dist(nrc*nsrc*kmax),&
        obst0(nrc*nsrc*kmax),periodRre(nrc*nsrc*kmax),stat=checkstat)
    obst=0
    obst0=0
    dist=0
    if(checkstat > 0)then
        write(6,*)'error with allocate'
    endif
    allocate(pvall(nrc*nsrc*kmax), stat=checkstat)
    if(checkstat > 0)then
        write(6,*)'error with allocate'
    endif

    istep=0
    istep2=0
    dall=0
    knumo=12345
    knum=0
    istep1=0
    pvall=0
    do
    read(87,'(a)',iostat=err) line
    if(err.eq.0) then
        if(line(1:1).eq.'#') then
            read(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
            if(wavetp.eq.2.and.veltp.eq.0) knum=period
            if(wavetp.eq.2.and.veltp.eq.1) stop 'can not deal with Rayleigh wave group data'
            if(wavetp.eq.1.and.veltp.eq.0) stop 'can not deal with Love wave phase data'
            if(wavetp.eq.1.and.veltp.eq.1) stop 'can not deal with Love wave group data'
            if(knum.ne.knumo) then
                istep=0
                istep2=istep2+1
            endif
            istep=istep+1
            istep1=0
            sta1_lat=(90.0-sta1_lat)*pi/180.0
            sta1_lon=sta1_lon*pi/180.0
            scxf(istep,knum)=sta1_lat
            sczf(istep,knum)=sta1_lon
            periods(istep,knum)=period
            wavetype(istep,knum)=wavetp
            igrt(istep,knum)=veltp
            nsrc1(knum)=istep
            knumo=knum
        else
            read(line,*) sta2_lat,sta2_lon,velvalue
            istep1=istep1+1
            dall=dall+1
            sta2_lat=(90.0-sta2_lat)*pi/180.0
            sta2_lon=sta2_lon*pi/180.0
            rcxf(istep1,istep,knum)=sta2_lat
            rczf(istep1,istep,knum)=sta2_lon
            call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
            dist(dall)=dist1
            periodRre(dall)=tRc(knum)
            obst0(dall)=dist1/velvalue
            pvall(dall)=velvalue
            nrc1(istep,knum)=istep1
        endif
    else
        exit
    endif
    enddo
    close(87)
    write(*,'(a,i7)') ' Number of all measurements',dall
    !--------------------------------------------------------------------------!
    !  Initialization
    allocate(depz(nz), stat=checkstat)
    maxvp = (nx-2)*(ny-2)*(nz-1)
    maxm =  (nx-2)*(ny-2)*(nz-1)*3
    maxnar = spfra*dall*nx*ny*nz*3 !sparsity fraction
    allocate(dv(maxm), stat=checkstat)
    allocate(norm(maxm), stat=checkstat)
    allocate(vsf(nx,ny,nz), stat=checkstat)
    allocate(vsftrue(nx,ny,nz), stat=checkstat)

    allocate(rw(maxnar), stat=checkstat)
    if(checkstat > 0)then
       write(6,*)'error with allocate: real rw'
    endif
    allocate(iw(2*maxnar+1), stat=checkstat)
    if(checkstat > 0)then
       write(6,*)'error with allocate: integer iw'
    endif
    allocate(col(maxnar), stat=checkstat)
    if(checkstat > 0)then
       write(6,*)'error with allocate:  integer iw'
    endif
    allocate(cbst(dall+maxm*maxm*3),dsyn(dall),stat=checkstat)
    ! allocate(deltaT(dall),stat=checkstat)
    allocate(sigmaT(dall),stat=checkstat)
    allocate(Tdata(dall),resbst(dall),stat=checkstat)

    allocate(Gctrue(nx-2,ny-2,nz-1),stat=checkstat)
    allocate(Gstrue(nx-2,ny-2,nz-1),stat=checkstat)
    allocate(obsTvs(dall),obsTaa(dall),stat=checkstat)
    allocate(Lsen_Gsc(nx*ny,kmaxRc,nz-1),STAT=checkstat)
    allocate(GVs(dall,maxvp),GGc(dall,maxvp),GGs(dall,maxvp),stat=checkstat)

    allocate(vsRela(nx-2,ny-2,nz-1),vsAbs(nx-2,ny-2,nz-1),stat=checkstat)
    allocate(vsReal(nx-2,ny-2,nz-1),stat=checkstat)
    allocate(gcf(nx-2,ny-2,nz-1),gsf(nx-2,ny-2,nz-1),stat=checkstat)
    allocate(iw2(2*maxnar+1), stat=checkstat)
    allocate( tRcV((nx-2)*(ny-2),kmaxRc),stat=checkstat)
    allocate(synT(dall),noise(dall),stat=checkstat)
    !--------------------------------------------------------------------------!
    !  READ INITIAL MODEL
    !--------------------------------------------------------------------------!
    !  Optional-Syn 1.1   MAKE SYNTHETIC TEST TO GENERATE OBSERVE DATA
    !--------------------------------------------------------------------------!
    write(*,*) ,'Forward Calculation Begin...'
    read(10,*)noiselevel
    write(*,'(a, f10.3)')'noise level: ',noiselevel
    write(66,'(a, f10.3)')'noise level: ',noiselevel

    inquire(file = 'MODVs.true', exist = ex)
    if (.not. ex) stop 'unable to open the MODVs.true'
    open(11,file='MODVs.true',status='old')
    vsftrue = 0
    read(11,*) (depz(i),i=1,nz)
    do k = 1,nz
        do j = 1,ny
            read(11,*) (vsftrue(i,j,k),i=1,nx)
        enddo
    enddo
    write(*,*) ' grid points in depth direction:(km)'
    write(*,'(50f6.2)') depz
    write(66,*) ' grid points in depth direction:(km)'
    write(66,'(50f6.2)') depz

    open(12,file='MODGc.true',status='old')
    open(13,file='MODGs.true',status='old')
    do k = 1,nz-1
        do j = 1,ny-2
            read(12,*) (Gctrue(i,j,k),i=1,nx-2)
            read(13,*) (Gstrue(i,j,k),i=1,nx-2)
        enddo
    enddo
    close(11)
    close(12)
    close(13)

    call CalRmax(nz,depz,minthk,rmax)
    !--------------------------------------------------------------------------!
    ! Optional-Syn 1.2: FwdObsTraveltime:  generate forward data using real model kernel and Gc,Gs model.
    ! forward traveltime based on the vsftrue.
    obsTvs=0
    obsTaa=0
    ! call FwdObsTraveltime(SenGTruefile, nx, ny, nz, maxvp, vsftrue,&
    ! Gctrue, Gstrue, obsTvs, obsTaa, dall, rmax, tRcV, Lsen_Gsc,&
    ! goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk,&
    ! scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, writepath)
    write(*,*)' Construct True Traveltime using Ture Sensitivity  Begin!'
    call FwdObsTraveltimeCPS(nx, ny, nz, maxvp, vsftrue,&
    Gctrue, Gstrue, obsTvs, obsTaa, dall, rmax, tRcV, Lsen_Gsc,&
    goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk,&
    scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, writepath)


    write(*,*)' Construct True Traveltime using True Sensitivity over!'
    ! forward calculate azimuthal  angle and amplitude
    open(42,file='period_Azm_tomo.real',status='replace',action='write')
    call FwdAzimuthalAniMap(nx,ny,nz,maxvp,&
        goxd,gozd,dvxd,dvzd,kmaxRc,tRc,&
        Gctrue,Gstrue,Lsen_Gsc,tRcV)
    ! add the gaussian noise only to anisotropic traveltime.
    sumObs=0.0
    synT=0
    noise=0
    sumNoise=0
    sumAdd=0
    do i=1,dall
        synT(i)=obsTvs(i)+obsTaa(i)
        ! noise(i)=synT(i)*gaussian()*noiselevel
        ! noiselevel need to be 0.5s
        noise(i)=gaussian()*noiselevel
        obst(i)=synT(i)+noise(i)
        sumObs=sumObs+abs(obsTaa(i))
        sumNoise=sumNoise+abs(noise(i))
        sumAdd=sumAdd+abs(noise(i)+obsTaa(i))
    enddo
    !--------------------------------------------------------------------------!
    !             Output the forward data
    !--------------------------------------------------------------------------!
    open(88,file='surfphase_forward.dat',action='write')
    count1=0
    sumTnos=0
    do knum=1,kmax
        do srcnum=1,nsrc1(knum)
            sta1_lat=scxf(srcnum,knum)
            sta1_lon=sczf(srcnum,knum)
            sta1_latD=90.0-sta1_lat*180.0/pi
            sta1_lonD=sta1_lon*180.0/pi
            write(88,'(a,2f11.6,3I3)')'#',sta1_latD,sta1_lonD,periods(srcnum,knum),wavetype(srcnum,knum),igrt(srcnum,knum)
            do istep=1,nrc1(srcnum,knum)
                sta2_lat=rcxf(istep,srcnum,knum)
                sta2_lon=rczf(istep,srcnum,knum)
                call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
                sta2_lat=90.0-sta2_lat*180.0/pi
                sta2_lon=sta2_lon*180.0/pi
                count1=count1+1
                Tvalue=obst(count1)
                velvalue=dist1/Tvalue

                velTrue=dist1/synT(count1)
                sumTnos=sumTnos+abs(velvalue-velTrue)/velTrue
                write(88,'(2f11.6,f9.5)')sta2_lat,sta2_lon,velvalue
            enddo
        enddo
    enddo
    close(88)

    open(88,file='Synthetic_fwd.dat')
    write(88,'(7a)'),'    Preiod        Distance(km)        T(s)       T_iso(s)    T_aa(s)    T_noe(s)     c(km/s)    c_iso(km/s)'
    do i=1,dall
        write(88,'(8f16.7)')periodRre(i),dist(i),obst(i),obsTvs(i),obsTaa(i),obst(i)-(obsTvs(i)+obsTaa(i)),dist(i)/obst(i),&
        dist(i)/obsTvs(i)
    enddo
    close(88)

    write(*,'(a,f13.3,a)')'  Max traveltime from Aniso: ',maxval(obsTaa(1:dall)),'s'
    ! write(*,*)'  location max obsTaa: ', MAXLOC(obsTaa(1:dall))
    write(*,'(a,f13.3,a)')'  Min traveltime from Aniso: ',minval(obsTaa(1:dall)),'s'
    write(*,'(a,f13.3,a)')'  Mean Abs t (s) from Aniso: ',sumObs/dall,'s'
    write(*,'(a,f13.3,a)')'  Mean Abs t (s) from Noise: ', sumNoise/dall,'s'
    write(*,'(a,f13.3,a)')'  Mean Abs t(Aniso+noisy) (s): ', sumAdd/dall,'s'
    write(*,'(a,f13.3)')'  Mean noisy Phase C (%): ', sumTnos/dall*100
    write(*,*)'--------------------make synthetic data over!-------------------------------'
    write(66,'(a,f13.3,a)')'  Max traveltime from Aniso: ',maxval(obsTaa(1:dall)),'s'
    write(66,'(a,f13.3,a)')'  Min traveltime from Aniso: ',minval(obsTaa(1:dall)),'s'
    write(66,'(a,f13.3,a)')'  Mean Abs t (s) from Aniso: ',sumObs/dall,'s'
    write(66,'(a,f13.3,a)')'  Mean Abs t (s) from Noise: ', sumNoise/dall,'s'
    write(66,'(a,f13.3,a)')'  Mean Abs t(Aniso+noisy) (s): ', sumAdd/dall,'s'
    write(66,'(a,f13.3)')  '  Mean noisy Phase C (%): ', sumTnos/dall*100

    !--------------------------------------------------------------------------!
    ! OUTPUT THE VELOCITY MODEL
    ! note: the input Vs model is grid model (nz point)
    ! but the inversion kernel and result is layer (nz-1 layer) model
    ! for Gc,Gs the input and result are both layer model, so the depth I set the lower depth as the output
    ! for Vs, input is grid model, while the inversion result is layered model. To compare two Vs models,
    ! I set the mid depth between the point as the depth index, and the real model use the average Vs.

    write(*,*)  'Program finishes successfully'
    write(66,*) 'Program finishes successfully'

    vsReal=0
    do k=1,nz-1
        do j=1,ny-2
            do i=1,nx-2
                vsref=(vsftrue(i+1,j+1,k)+vsftrue(i+1,j+1,k+1))/2
                vsReal(i,j,k)=vsref
            enddo
        enddo
    enddo
    open(71,file='Gc_Gs_model.real')
    call WriteAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,Gctrue,Gstrue,vsReal,71)
    close(71)

    open(72,file='Vs_model.real')
    do k=1,nz-1
        do j=1,ny-2
            do i=1,nx-2
                vsref=(vsftrue(i+1,j+1,k)+vsftrue(i+1,j+1,k+1))/2
                write(72,'(5f10.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,(depz(k)+depz(k+1))/2,vsref
            enddo
        enddo
    enddo
    close(72)


    endT=OMP_get_wtime()
    write(*,'(a,f13.0,a)')"     All time cost=",endT-startT,"s"

    write(*,*),'Output True velocity model to Vs_model.real'

    write(*,*),'Output inverted shear velocity model &
          to Vs_model_Syn.rela and Vs_model_Syn.abs'
    write(66,*),'Output True Gc Gs model &
        to Gc_model.real Gs_model.real'
    write(66,*),'Output inverted shear velocity model &
        to Vs_model_Syn.rela and Vs_model_Syn.abs'

    close(10) ! surfaniso.in
    close(66) !close surf_tomo.log

    deallocate(obst)
    deallocate(dsyn)
    deallocate(dist,periodRre)
    deallocate(pvall)
    deallocate(depz)
    deallocate(scxf,sczf)
    deallocate(rcxf,rczf)
    deallocate(wavetype,igrt,nrc1)
    deallocate(nsrc1,periods)
    deallocate(rw)
    deallocate(iw,col)
    deallocate(cbst)
    deallocate(dv)
    deallocate(norm)
    deallocate(vsf)
    deallocate(vsftrue)
    deallocate(sigmaT,Tdata,resbst)
    deallocate(Gctrue,Gstrue)
    deallocate(obsTvs,obsTaa)
    deallocate(GVs,GGs,GGc,Lsen_Gsc)
    deallocate(vsRela,vsAbs)
    deallocate(vsReal)
    deallocate(gcf,gsf)
    deallocate(iw2)
    if(kmaxRc.gt.0) deallocate(tRc)
    deallocate(tRcV)
    deallocate(synT,noise)
end program



subroutine WriteAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,Gc,Gs,vsBg,Idout)
    integer nx,ny,nz
    real goxd,gozd
    real dvxd,dvzd
    real depz(nz)
    real Gc(nx-2,ny-2,nz-1), Gs(nx-2,ny-2,nz-1), vsBg(nx-2,ny-2,nz-1)
    integer Idout
    integer:: k,j,i
    real cosTmp,sinTmp
    real*8 :: pi=3.1415926535898
    do k=1,nz-1
        do j=1,ny-2
            do i=1,nx-2
                cosTmp=Gc(i,j,k)
                sinTmp=Gs(i,j,k)
                AzimAmp=0.5*sqrt(cosTmp**2+sinTmp**2)
                AzimAng=atan2(sinTmp,cosTmp)/pi*180
                if (AzimAng.LT.0.0) AzimAng=AzimAng+360
                AzimAng=0.5*AzimAng
                ! format: lon lat  depth(lower) Vs(mid-depth)  Angle  Amp(%/100) Gc/L(%)  Gs/L(%)
                write(Idout,'(8f10.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k+1),vsBg(i,j,k),AzimAng,AzimAmp,&
                   Gc(i,j,k)*100,Gs(i,j,k)*100
            enddo
        enddo
    enddo
end subroutine


subroutine CalRmax(nz,depz,minthk0,rmax)
    implicit none
    integer nz
    real depz(nz)
    real minthk0,minthk
    integer rmax
    integer k,i
    integer NL
    parameter (NL=200)
    integer nsublay(NL)
    real thk
    rmax=0
    do i=1,nz-1
        thk = depz(i+1)-depz(i)
        minthk = thk/minthk0
        nsublay(i) = int((thk+1.0e-4)/minthk) + 1
    enddo
    rmax=sum(nsublay(1:nz-1))
end subroutine
