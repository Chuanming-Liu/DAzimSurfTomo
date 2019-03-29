      ! CODE FOR SURFACE WAVE TOMOGRAPHY USING DISPERSION MEASUREMENTS
      ! VERSION:
      !      1.0
        program SurfAniso
        use lsmrModule, only:lsmr
	use  lsmrblasInterface, only : dnrm2
	use omp_lib
        implicit none

! VARIABLE DEFINE
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
            real weightVs,weight0,weight1,weightGcs,weight
            real damp,dampAA,dampVs
            real minthk
            integer kmax,kmaxRc
            real*8,dimension(:),allocatable:: tRc
            real,dimension(:),allocatable:: depz
            integer itn
            integer nout
            integer localSize
            real mean,std_devs,balances,balanceb
            integer msurf
            real,dimension(:),allocatable:: obst,dsyn,cbst,dist
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
                real,dimension(:),allocatable::dv
                real,dimension(:,:,:),allocatable::vsf
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
            real acond,anorm,arnorm,rnorm,xnorm
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
            integer writepath
            real averdws
            real maxnorm
! Added by C.M. L
            character SenGTruefile*200, SenGfile*200
            integer lay
            integer maxm
            real, dimension(:), allocatable:: sigmaT,resSigma
             real, dimension(:), allocatable:: Tdata,resbst,fwdT,fwdTvs,fwdTaa,Tref,RefTaa,resbst_iso

            real meandeltaT,stddeltaT
            real guassVs,LcorrXY
            real sigmaDeltaVs,sigmaGc,sigmaGs
            real sigmaGcs
             real LcorrZZ
             real,dimension(:,:),allocatable:: GVs,GGs,GGc
             real,dimension(:,:,:), allocatable :: Lsen_Gsc
            real gaussian
            external gaussian
            real synT
!            real,dimension(:,:,:),allocatable::vsRela,vsAbs,vsReal
            real,dimension(:,:,:),allocatable:: gcf,gsf
            integer nar1,nar2
            real vsref
            REAL res2Norm,Mnorm2,Enorm2,Vsnorm2,Gcnorm2,Gsnorm2
            REAL resNSigma2Norm
            REAL Enorm2Lame,Vsnorm2Lame,Gcnorm2Lame,Gsnorm2Lame
            REAL MwNorm2
            INTEGER writeperiod
            REAL*8,dimension(:,:),allocatable:: tRcV
            REAL sumObs
            INTEGER Refwritepath,Realwritepath
            REAL VariGc,VariGs,VariVs
            INTEGER count4,para
            REAL LambdaVs,LambdaGc,LambdaGs
            REAL PreRes,meanAbs
            INTEGER iter_mod,iso_inv,rmax
            CHARACTER(LEN=2) :: id='00'
            CHARACTER(LEN=60) :: filename
            REAL*8 :: startT,endT
            LOGICAL:: isTikh
! OPEN FILES FIRST TO OUTPUT THE PROCESS
            startT=OMP_get_wtime()
            nout=36
            OPEN(nout,file='lsmr.txt')
! OUTPUT PROGRAM INFOMATION
            WRITE(*,*)
            WRITE(*,*),'                         SurfAniso'
            WRITE(*,*)'PLEASE contact Chuanming Liu &
                (chuanmingliu@foxmail.com) IF you find any bug.'
            WRITE(*,*)

! READ INPUT FILE
            IF (iargc() < 1) THEN
              WRITE(*,*) 'input file [SurfAniso.in(default)]:'
              READ(*,'(a)') inputfile
              IF (len_trim(inputfile) <=1 ) THEN
                  inputfile = 'SurfAniso.in'
              ELSE
                  inputfile = inputfile(1:len_trim(inputfile))
              ENDIF
            ELSE
                CALL getarg(1,inputfile)
            ENDIF
            inquire(file = inputfile, exist = ex)
            IF (.not. ex)   STOP 'unable to OPEN the inputfile'

            OPEN(10,file=inputfile,status='old',action='READ')
            READ(10,'(a30)')dummy
            READ(10,'(a30)')dummy
            READ(10,'(a30)')dummy
            READ(10,*)SenGfile
            READ(10,*)datafile
            READ(10,*) nx,ny,nz
            READ(10,*) goxd,gozd
            READ(10,*) dvxd,dvzd
            READ(10,*) nsrc
            READ(10,*) minthk
            READ(10,*) maxiter
            READ(10,*) spfra
            READ(10,'(a30)')dummy
            READ(10,*)weight1
            READ(10,*)weight0
             READ(10,*)dampVs
            READ(10,*)dampAA
            READ(10,*) isTikh
            READ(10,*) LcorrXY,LcorrZZ

            WRITE(*,*) 'integrated  sensitivity kernel for Vs and Gcs file based on MOD  '
            WRITE(*,'(a)')SenGfile
            WRITE(*,*)'input Rayleigh wave phase velocity data file:'
            WRITE(*,'(a)')datafile
            WRITE(*,*)  'model origin:latitude,longitue'
            WRITE(*,'(2f10.4)') goxd,gozd
            WRITE(*,*) 'grid spacing:latitude,longitue'
            WRITE(*,'(2f10.4)') dvxd,dvzd
            WRITE(*,*) 'model dimension:nx,ny,nz'
            WRITE(*,'(3i5)') nx,ny,nz
            WRITE(*,*)'depth refined interval layer '
            WRITE(*,'(f8.1)')minthk
            WRITE(*,*)'weight for Vs '
            WRITE(*,'(f8.1)')weight1
            WRITE(*,*)'weight for Gc, Gs '
            WRITE(*,'(f8.1)')weight0
            WRITE(*,*)'damp for Vs '
            WRITE(*,'(f8.1)')dampVs
            WRITE(*,*)'damp for Gc, Gs '
            WRITE(*,'(f8.1)')dampAA
            WRITE(*,*)' Regularization Type: (T) 1st order Tikhonov ;(F) Gaussian'
            ! WRITE(*,*)isTikh



            WRITE(*,*) 'correlation length: XY, ZZ (km)'
            WRITE(*,'(f8.1)')LcorrXY, LcorrZZ

            IF (nz.LE.1)  STOP 'error nz value.'

            READ(10,*) sigmaDeltaVs
            WRITE(*,*)'model perturbation: Vs (km/s)'
            WRITE(*,'(50f9.5)')sigmaDeltaVs

            READ(10,*) sigmaGcs
            WRITE(*,*)'model perturbation: Gc,Gs (%) layer'
            WRITE(*,'(50f9.5)')sigmaGcs*100

            READ(10,'(a30)')dummy
            READ(10,*) kmaxRc
            WRITE(*,*) 'number of period'
            WRITE(*,'(i6)') kmaxRc

            IF(kmaxRc.gt.0)THEN
               ALLOCATE(tRc(kmaxRc),STAT=checkstat)
               IF (checkstat > 0) STOP 'error allocating RP'
               READ(10,*)(tRc(i),i=1,kmaxRc)
            ELSE
                STOP 'Can only deal with Rayleigh wave phase velocity data!'
            ENDIF

	    WRITE(logfile,'(a,a)')trim(inputfile),'.log'
!            OPEN(66,file=logfile,action='WRITE')
            OPEN(66,file=logfile)
            WRITE(66,*)
            WRITE(66,*),'                    SurfAniso'
            WRITE(66,*)'PLEASE contact Chuanming Liu &
                 (chuanmingliu@foxmail.com) IF you find any bug.'
            WRITE(66,*)
            WRITE(66,*) 'model origin:latitude,longitue'
            WRITE(66,'(2f10.4)') goxd,gozd
            WRITE(66,*) 'grid spacing:latitude,longitue'
            WRITE(66,'(2f10.4)') dvxd,dvzd
            WRITE(66,*) 'model dimension:nx,ny,nz'
            WRITE(66,'(3i5)') nx,ny,nz

            WRITE(*,*)'Rayleigh wave phase velocity used,periods:(s)'
            WRITE(*,'(50f6.2)')(tRc(i),i=1,kmaxRc)
            WRITE(66,*)'Rayleigh wave phase velocity used,periods:(s)'
            WRITE(66,'(50f6.2)')(tRc(i),i=1,kmaxRc)
            WRITE(66,'(a)')' --------------------%%%%%%%%%%%%%%%%%%%----------------------------'

!            CLOSE(10)
            nrc=nsrc
            kmax=kmaxRc
            lay=nz-1
!------------------------------------------------------------------------------------------------!
! READ MEASUREMENTS
            inquire(file = datafile, exist = ex)
            IF (.not. ex) THEN
              WRITE(66,'(a)')'unable to OPEN the datafile'
              CLOSE(66)
              STOP 'unable to OPEN the datafile'
             ENDIF
             WRITE(*,*)'begin load data file.....'

            OPEN(unit=87,file=datafile,status='old')
            ALLOCATE(scxf(nsrc,kmax),sczf(nsrc,kmax),&
            rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),STAT=checkstat)
            IF(checkstat > 0)  WRITE(6,*)'scxf error with ALLOCATE','nsrc=',nsrc,'kmax=',kmax

            ALLOCATE(periods(nsrc,kmax),wavetype(nsrc,kmax),nrc1(nsrc,kmax),nsrc1(kmax),&
            igrt(nsrc,kmax),STAT=checkstat)
            IF(checkstat > 0) WRITE(6,*)'error with ALLOCATE'

            ALLOCATE(obst(nrc*nsrc*kmax),dist(nrc*nsrc*kmax),STAT=checkstat)
            obst=0
            dist=0
            IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE'
           ALLOCATE(pvall(nrc*nsrc*kmax),STAT=checkstat)
            IF(checkstat > 0)    WRITE(6,*)'error with ALLOCATE'

            istep=0
            istep2=0
            dall=0
            knumo=12345
	    knum=0
	    istep1=0
	    pvall=0
            DO
               READ(87,'(a)',iostat=err) line
               IF(err.eq.0) THEN
                       IF(line(1:1).eq.'#') THEN
                          READ(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
                          IF(wavetp.eq.2.and.veltp.eq.0) knum=period
                          IF(wavetp.eq.2.and.veltp.eq.1) STOP 'can not deal with Rayleigh wave group data'
                          IF(wavetp.eq.1.and.veltp.eq.0) STOP 'can not deal with Love wave phase data'
                          IF(wavetp.eq.1.and.veltp.eq.1) STOP 'can not deal with Love wave group data'
                          IF(knum.ne.knumo) THEN
                             istep=0
                             istep2=istep2+1
                          ENDIF
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
                      ELSE
                          READ(line,*) sta2_lat,sta2_lon,velvalue
                          istep1=istep1+1
                          dall=dall+1
                          sta2_lat=(90.0-sta2_lat)*pi/180.0
                          sta2_lon=sta2_lon*pi/180.0
                          rcxf(istep1,istep,knum)=sta2_lat
                          rczf(istep1,istep,knum)=sta2_lon
                         CALL delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
                         dist(dall)=dist1
                         obst(dall)=dist1/velvalue
                         pvall(dall)=velvalue
                         nrc1(istep,knum)=istep1
                      ENDIF
               ELSE
                  EXIT
               ENDIF
            ENDDO
            CLOSE(87)
            WRITE(*,'(a,i7)') ' Number of all measurements',dall
!---------------------------------------------------------------------------!
!  Initialization
            ALLOCATE(depz(nz), STAT=checkstat)
            ALLOCATE(vsf(nx,ny,nz), STAT=checkstat)

            maxvp = (nx-2)*(ny-2)*(nz-1)
            maxm =  (nx-2)*(ny-2)*(nz-1)*2
            maxnar = spfra*dall*nx*ny*nz*2 !sparsity fraction

            ALLOCATE(dv(maxm), stat=checkstat)
      	    ALLOCATE(rw(maxnar), STAT=checkstat)
            ALLOCATE(iw(2*maxnar+1), STAT=checkstat)
            ALLOCATE(col(maxnar), STAT=checkstat)
            !ALLOCATE(cbst(dall+maxm*maxm),dsyn(dall),STAT=checkstat)

            ALLOCATE(cbst(dall+maxm*2),dsyn(dall),STAT=checkstat)


            ALLOCATE(dsyn(dall),STAT=checkstat)
            ALLOCATE(sigmaT(dall),Tref(dall),STAT=checkstat)
            IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE:  sigma,Tref'


            ALLOCATE(Tdata(dall),resbst(dall),fwdT(dall),fwdTvs(dall), fwdTaa(dall),RefTaa(dall),STAT=checkstat)
             ALLOCATE(resbst_iso(dall),STAT=checkstat)
             IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE:  Tdata,resbst'
            ALLOCATE(resSigma(dall),STAT=checkstat)
            ALLOCATE(Lsen_Gsc(nx*ny,kmaxRc,nz-1),STAT=checkstat)
            ALLOCATE(GVs(dall,maxvp),GGc(dall,maxvp),GGs(dall,maxvp),STAT=checkstat)
            IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE:  GVs,GGc,GGs'
            ALLOCATE(gcf(nx-2,ny-2,nz-1),gsf(nx-2,ny-2,nz-1),STAT=checkstat)
             IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE:  gcf,gsf'
            ALLOCATE( tRcV((nx-2)*(ny-2),kmaxRc),STAT=checkstat)
            IF(checkstat > 0)  WRITE(6,*)'error with ALLOCATE:  tRcV'
!-----------------------------------------------------------------------------------------------------!
!  READ INITIAL MODEL
            OPEN(11,file='MOD',status='old')
            vsf=0
            READ(11,*) (depz(i),i=1,nz)
            DO k = 1,nz
                DO j = 1,ny
                   READ(11,*)(vsf(i,j,k),i=1,nx)
                ENDDO
            ENDDO
            CLOSE(11)
            WRITE(*,*) ' grid points in depth direction:(km)'
            WRITE(*,'(50f6.2)') depz

            CALL  CalRmax(nz,depz,minthk,rmax)
!-----------------------------------------------------------------------------------------------------!
!                     Iteration Part  ITERATE UNTILL CONVERGE
!-----------------------------------------------------------------------------------------------------!
!  IF iter is  even: invert for isotropic model Vs
!  IF iter is odd: invert for  anisotropic model Gc/L,Gs/L
!-----------------------------------------------------------------------------------------------------!
            WRITE(6,*)'----------------------INVERSION BEGIN--------------------------------'
            WRITE(66,*)'----------------------INVERSION BEGIN--------------------------------'
            RefTaa=0
            DO iter = 1,maxiter

               iter_mod=mod(iter,2)
            WRITE(6,*)'!-------------%%%%%%%%%%%%%%%%---------------------------!'
            WRITE(66,*)'!-------------%%%%%%%%%%%%%%%%---------------------------!'
               IF (iter_mod.eq.0) THEN
                   iso_inv=1
                   WRITE(66,*)iter,'th iteration, invert for isotropic Vs para.'
                   WRITE(6,*)iter,'th iteration, invert for isotropic Vs para.'
                    maxm =  maxvp
               ELSE
                    iso_inv=0
                    WRITE(66,*)iter,'th iteration, invert for anisotropic Gc, Gs para.'
                    WRITE(6,*)iter,'th iteration, invert for anisotropic Gc, Gs para.'
                    maxm =  maxvp*2
               ENDIF

! COMPUTE SENSITIVITY MATRIX
               WRITE(66,'(a)')' --------------------%%%%%%%%%%%%%%%%----------------------------'
               WRITE(6,'(a)')' --------------------%%%%%%%%%%%%%%%%----------------------------'
               Refwritepath = 0 ! switch of WRITE ray path
               WRITE(*,*) iter,'th iteration...',' computing sensitivity matrix.'
               WRITE(66,*)iter,'th iteration...','computing sensitivity matrix'
               dsyn=0
               GGc=0
               GGs=0
               GVs=0
               tRcV=0
               iw = 0
               rw = 0.0
               col = 0
               IF (iso_inv.eq.0)THEN
                  CALL CalSurfGAniso(SenGfile,nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
                   GGc,GGs,Lsen_Gsc,dall,rmax,tRcV,&
                   goxd,gozd,dvxd,dvzd,kmaxRc,tRc,periods,depz,minthk,&
                   scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax, nsrc,nrc,nar,Refwritepath)
               ELSEIF(iso_inv.eq.1)THEN
                   CALL CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
                    GVs,dall,&
                    goxd,gozd,dvxd,dvzd,kmaxRc,tRc,periods,depz,minthk,&
                     scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,nsrc,nrc,nar)
               ENDIF
              WRITE(*,*) ' Finish G matrix calculation.'
!-----------------------------------------------------------------------------------------------------!
! output the MOD corresponding phase velocity map.
               IF (iter.eq.1)THEN
                OPEN(77,file='period_phaseVMOD.dat')
                WRITE(66,*)'output the MOD corresponding period Tomo map-period_phaseVMOD.dat'
                CALL WTPeriodPhaseV(nx,ny,gozd,goxd,dvzd,dvxd,kmaxRc,tRc,tRcV,77)
              END IF
!-----------------------------------------------------------------------------------------------------!
! CALCULATE DATA RESIDUAL (Delta T)
! cbst: travel time residual.
! deltaT: relative travel time residual
! sigmaT: Cd(i,i)=sigmaT**2
! opt: (Cd^-1)^0.5-->1/sigmaT

! dsyn: predicated isotropic travel time
               WRITE(*,*),iter,'th iteration...',' computing reference model data residual.'
               WRITE(66,*),iter,'th iteration...',' computing reference model data residual.'
               cbst=0
               Tdata=0
               sigmaT=0
               Tref=0
               IF (iso_inv.eq.0)THEN
                   DO i=1,dall
                      Tref(i)=dsyn(i)
	                    cbst(i) = obst(i) - Tref(i)
                      Tdata(i)=cbst(i)
                    ENDDO
               ELSEIF(iso_inv.eq.1)THEN
                    DO i = 1,dall
   	              Tref(i)=dsyn(i)+RefTaa(i)
   	              cbst(i) = obst(i) - Tref(i)
                      Tdata(i)=cbst(i)
	           ENDDO
               ENDIF

! Statistic of the data residual using the reference model.
               mean = sum(cbst(1:dall))/dall
               std_devs = sqrt(sum((cbst(1:dall)-mean)**2)/dall)
               meanAbs=sum(abs(cbst(1:dall)))/dall

               WRITE(*,'(i2,a)')  ,iter,'th Inversion.--Input Data Space'

               WRITE(66,'(i2,a)')  ,iter,'th Inversion.--Input Data Space'
               IF (iso_inv.eq.0)THEN
                   WRITE(66,'(a,f12.4,a,f12.4,a,f12.4)'),' Whole  Traveltime Residual (Tobs-Tiso) mean and std_devs of &
                   residual: ',mean*1000,'ms ',1000*std_devs,'ms '
                   WRITE(66,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
                   WRITE(6,'(a,f12.4,a,f12.4,a,f12.4)'),' Whole  Traveltime Residual (Tobs-Tiso) mean and std_devs of &
                   residual: ',mean*1000,'ms ',1000*std_devs,'ms '
                   WRITE(6,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
               ELSEIF(iso_inv.eq.1)THEN
                  WRITE(66,'(a,f12.4,a,f12.4,a,f12.4)'),' Whole  Traveltime Residual (Tobs-Tiso-Taa) mean and std_devs of &
                   residual: ',mean*1000,'ms ',1000*std_devs,'ms '
                   WRITE(66,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
                  WRITE(6,'(a,f12.4,a,f12.4,a,f12.4)'),' Whole  Traveltime Residual (Tobs-Tiso-Taa) mean and std_devs of &
                   residual: ',mean*1000,'ms ',1000*std_devs,'ms '
                   WRITE(6,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
               ENDIF
 !-----------------------------------------------------------------------------------------------------!
! Calculate the sigma for the constraint on the data used in the inversion.
	       CALL CalDdatSigma(dall,obst,cbst,sigmaT,meandeltaT,stddeltaT)
	       DO i=1,dall
	            cbst(i)=cbst(i)*1/sigmaT(i)
               ENDDO
               WRITE(6,*),iter,'th iteration...','Mean sigma',sum(sigmaT(1:dall))/dall
               WRITE(66,*),iter,'th iteration...','Mean sigma',sum(sigmaT(1:dall))/dall
               WRITE(66,*),iter,'th iteration...',' Mean relative abs traveltime residual(dt/t0) %',meandeltaT*100
               WRITE(66,*),iter,'th iteration...',' Std relative traveltime residual(dt/t0)',stddeltaT
                WRITE(6,*),iter,'th iteration...',' Mean relative abs traveltime residual(dt/t0) %',meandeltaT*100
               WRITE(6,*),iter,'th iteration...',' Std relative traveltime residual(dt/t0)',stddeltaT
! ADDING THE Cd^-1 TO G
	       DO i = 1,nar
		         rw(i) = rw(i)*1/sigmaT(iw(1+i))
               ENDDO
! WRITE OUT RESIDUAL FOR THE FIRST AND LAST ITERATION
               IF(iter.eq.1) THEN
                  OPEN(88,FILE='Traveltime_use.dat')
                  WRITE(88,'(7a)'),'Dist(km)        T_obs(s)       T_forward(s)         Res(s)'
                  DO i=1,dall
                      WRITE(88,*) dist(i),obst(i),dsyn(i),Tdata(i)
                  ENDDO
                  CLOSE(88)
               ENDIF
!-----------------------------------------------------------------------------------------------------!
! Setting with iteration.
! After 1st iteration, use the res2Norm and Mnorm2 to set weight value in the next iteration.
!-----------------------------------------------------------------------------------------------------!
!  ADDING GUASSIAN REGULARIZATION
  !  A*x = b
  ! m       input      m, the number of rows in A.
  ! n       input      n, the number of columns in A.'
  ! count3: data, d, counter index
  ! nar: number of G which value is not zero.IF Gn*m does not contain zeros,  nar =n*m
  ! rw: G none zero value
  ! col: counter index on model parameter, col(i) means model parameter order at i th data in G ---m index
  ! iw: counter index on data.---n index, iw(i) means data order at i th data in G, i=1: nar
!-----------------------------------------------------------------------------------------------------!
               WRITE(*,*)iter,'th iteration...',' computing guassian regularization ...'
               nar1=nar
!               weightGcs=dnrm2(dall,cbst,1)/dall*weight0
!               weightVs=dnrm2(dall,cbst,1)/dall*weight1
               weightGcs=weight0
               weightVs=weight1
               count3=0

               IF ( isTikh )THEN
!                  weightVs=dnrm2(dall,cbst,1)/dall*weight1
                   CALL TikhonovRegularization(nx,ny,nz,maxvp,dall,nar,rw,iw,col,count3,iso_inv,weightGcs,weightVs)
               ELSE
                   CALL GaussianLS(nz,ny,nx,nar,maxvp,dall ,gozd,goxd,dvzd,dvxd,depz,LcorrXY,LcorrZZ,&
                               weightVs,weightGcs,sigmaDeltaVs,sigmaGcs,rw,iw,col,count3,iso_inv)
               ENDIF

               FORALL(i=1:count3)
                   cbst(dall+i)=0
               END FORALL

               !---------------------------------------------------------@


               m = dall + count3
               n=maxm


               iw(1)=nar
               DO i=1,nar
                   iw(1+nar+i)=col(i)
               ENDDO

               iF (nar > maxnar) STOP 'increase sparsity fraction(spfra)'
               WRITE(*,'(i2,a)'),iter,'th iteration...'
               WRITE(*,'(a,i9)') '  Model number=',n
               WRITE(*,'(a,i9)') '  Model Space maxvp=',maxvp
               WRITE(*,'(a,i9)') '  Data number=',dall

               WRITE(*,'(a,i9,a)')'The matrix  A  has',m,'rows'
               WRITE(*,'(a,i9)') '  Regularization matrix row number=',count3
               WRITE(*,'(a,i9)')'  nar (G none zeros)=', nar
               WRITE(*,'(a,i9)') '  L row number=',(nar-nar1)/2
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
! CALLING IRLS TO SOLVE THE PROBLEM
               WRITE(*,*)'-----------------------------------------------------------'
               WRITE(*,*)' solving the problem by LSMR ...'
               leniw = 2*nar+1
               lenrw = nar
               dv = 0
               atol = 1e-6
               btol = 1e-6
               conlim = 2000
               itnlim = 1000
               istop = 0
               anorm = 0.0
               acond = 0.0
               arnorm = 0.0
               xnorm = 0.0
               localSize = n/4



               WRITE(*,'(a,f7.1)')'  LcorrXY in L (km):',LcorrXY
               WRITE(*,'(a,f7.1)')'  LcorrZZ in L (km):',LcorrZZ


               WRITE(66,'(a,f7.1)')'  LcorrXY in L (km)=',LcorrXY
               WRITE(66,'(a,f7.1)')'  LcorrZZ in L (km)=',LcorrZZ

               IF (iso_inv.eq.0)THEN
                  damp=dampAA
                   WRITE(*,'(a,f8.3)')'  damp for Gcs:',dampAA
                   WRITE(66,'(a,f8.1)')'  damp for Gcs=',dampAA
                   WRITE(*,'(a,2f8.3)')' Gcs weight and weight0 for L=',weightGcs,weight0
                   WRITE(66,'(a,2f8.3)')' Gcs weight and weight0 for L=',weightGcs,weight0
                    WRITE(*,'(a,20f8.5)')'  sigma Gc, Gs in L(%):',sigmaGcs*100
                    WRITE(66,'(a,20f8.4)')'  sigma Gc, Gs in L (%):',sigmaGcs*100
               ELSEIF (iso_inv.eq.1)THEN
                   damp=dampVs
                   WRITE(*,'(a,f8.3)')'  damp for dVs:',dampVs
                   WRITE(66,'(a,f8.1)')'  damp for dVs=',dampVs
                   WRITE(*,'(a,2f8.3)')' Vs weight and weight0 for L=',weightVs,weight1
                   WRITE(66,'(a,2f8.3)')'  Vs weight and weight0 for L=',weightVs,weight1
                    WRITE(*,'(a,20f8.5)')'  sigmaVs in L (km/s):',sigmaDeltaVs
                    WRITE(66,'(a,20f8.4)')'  sigmaVs in L(km/s):',sigmaDeltaVs
               ENDIF

               CALL LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
               atol, btol, conlim, itnlim, localSize, nout,&
               dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
               IF(istop==3) print*,'istop = 3, large condition number'
               WRITE(*,*)'Finish LSMR.......'
               WRITE(*,*)'istop=',istop
               WRITE(*,*)'anorm=',anorm
               WRITE(*,*)'acond=',acond
               WRITE(*,*)'rnorm=',rnorm
               WRITE(*,*)'arnorm=',arnorm
               WRITE(*,*)'xnorm=',xnorm

!-----------------------------------------------------------------------------------------------------!
! CONSTRUCT THE VELOCITY MODEL
! For vsRela will be used in subroutine: CalSurfGAniso, it will be important to get a right value.
! vsf will not change, which is from MOD.
               IF (iso_inv.eq.0)THEN
                   gcf=0
                   gsf=0
                   DO k=1,nz-1
                       DO j=1,ny-2
                          DO i=1,nx-2
                              gcf(i,j,k)=dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                              gsf(i,j,k)=dv(maxvp+(k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                          ENDDO
                       ENDDO
                   ENDDO
               ELSEIF(iso_inv.eq.1)THEN
                   DO k=1,nz-1
                       DO j=1,ny-2
                          DO i=1,nx-2
                             vsf(i+1,j+1,k)=vsf(i+1,j+1,k)+dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)

                          ENDDO
                       ENDDO
                   ENDDO
               ENDIF
!-----------------------------------------------------------------------------------------------------!
! INVERSION RESULT
               IF (iso_inv.eq.0)THEN
                    WRITE(6,'(a,3f10.4)'),'  min  max and abs mean  Gc/L (%) ',&
                    minval(dv(1:maxvp))*100,maxval(dv(1:maxvp))*100,&
                    sum(abs(dv(1:maxvp)))/maxvp*100
                    WRITE(6,'(a,3f10.4)'),'   min  max and abs mean  Gs/L (%) ',&
                    minval(dv(maxvp+1:maxvp+maxvp))*100,maxval(dv(maxvp+1:maxvp+maxvp))*100,&
                    sum(abs(dv(maxvp+1:maxvp+maxvp)))/maxvp*100
                   WRITE(66,'(a,3f10.4)'),'  min  max and abs mean Gc/L (%)   ',&
                    minval(dv(1:maxvp))*100,maxval(dv(1:maxvp))*100,&
                    sum(abs(dv(1:maxvp)))/maxvp*100
                    WRITE(66,'(a,3f10.4)'),'  min  max and abs mean Gs/L (%)   ',&
                    minval(dv(maxvp+1:maxvp+maxvp))*100,maxval(dv(maxvp+1:maxvp+maxvp))*100,&
                    sum(abs(dv(maxvp+1:maxvp+maxvp)))/maxvp*100
                   WRITE(66,'(a)'),''
                    DO k=1,nz-1
                        VariGc=sum(abs(gcf(1:nx-2,1:ny-2,k)))/((nx-2)*(ny-2))
                        VariGs=sum(abs(gsf(1:nx-2,1:ny-2,k)))/((nx-2)*(ny-2))
                         WRITE(66,'(a,f4.1,a,f4.1,a,2f10.3)') '  Depth= ',depz(k),'--',depz(k+1),'km  Abs Mean Variation (%)  Gc &
                         Gs',VariGc*100,VariGs*100
                          WRITE(6,'(a,f4.1,a,f4.1,a,2f10.3)') '  Depth= ',depz(k),'--',depz(k+1),'km  Abs Mean Variation (%)  Gc &
                         Gs',VariGc*100,VariGs*100
                    ENDDO
               ELSEIF(iso_inv.eq.1)THEN
                     WRITE(6,'(a,3f10.4)'),'  min  max and abs mean  dVs(km/s)',&
                     minval(dv(1:maxvp)),maxval(dv(1:maxvp)),sum(abs(dv(1:maxvp)))/maxvp
                     WRITE(66,'(a,3f10.4)'),'  min  max and abs mean dVs(km/s) ',&
                     minval(dv(1:maxvp)),maxval(dv(1:maxvp)),sum(abs(dv(1:maxvp)))/maxvp
                     WRITE(66,'(a)'),''
                    DO k=1,nz-1
                        VariVs=sum(abs(dv((k-1)*(nx-2)*(ny-2)+1:k*(nx-2)*(ny-2) )))/((nx-2)*(ny-2))
                        WRITE(66,'(a,f4.1,a,f4.1,a,f10.3)') '  Depth= ',depz(k),'--',depz(k+1),&
                        'km  Abs Mean Variation (km/s) Vs',VariVs
                        WRITE(6,'(a,f4.1,a,f4.1,a,f10.3)') '  Depth= ',depz(k),'--',depz(k+1),&
                        'km  Abs Mean Variation (km/s) Vs',VariVs
                    ENDDO
               ENDIF
!-----------------------------------------------------------------------------------------------------!
! Calculate ||Lm||2
! same as GaussianLS
               Mnorm2=0
               MwNorm2=0
               IF (iso_inv.eq.0)THEN
                  weight=weightGcs
               ELSEIF(iso_inv.eq.1)THEN
                 weight=weightVs
               ENDIF
               IF (weight.NE.0)THEN
                  count3=nar-nar1
                  CALL Calmodel2Norm(nar1,nar,maxvp,count3,rw,col,dv,weight,Mnorm2,MwNorm2,isTikh)
                  WRITE(*,*)'Calmodel2Norm finished!'
              ENDIF
!-----------------------------------------------------------------------------------------------------!
! Calculate ||W(Gm-d)||2
! fist forward calculate the traveltime.
! Tdata: inversion used data: tobs- tref(iso)
! Calculate deltaT=deltaTvs+deltaTaa
! Here restT=deltaT_in-deltaT_out---use the old GGc GGs
               resbst=0
               res2Norm=0
               resNSigma2Norm=0
               WRITE(id,'(I2.2)') iter
               filename='Traveltime_Delta_'//TRIM(id)//'th.dat'
               OPEN(88,file=TRIM(filename))

               IF (iso_inv.eq.0)THEN
                   fwdTaa=0
                   CALL  CalGcsReslNorm(maxvp,dall,GGc,GGs,dv,sigmaT,Tdata,fwdTaa,resbst,res2Norm,resNSigma2Norm,PreRes)
                    WRITE(88,'(a)'),' Distance(km)        T_obs(s)         T_ref-iso      Data_Res&
                         Delta_aa(s)              Res(s)'
                    DO i=1,dall
                        WRITE(88,*) dist(i),obst(i),dsyn(i),Tdata(i),fwdTaa(i),resbst(i)
                    ENDDO
                    CLOSE(88)
               ELSEIF(iso_inv.eq.1)THEN
                    CALL CalVsReslNorm(maxvp,dall,GVs,dv,sigmaT,Tdata,fwdTvs,resbst,res2Norm,resNSigma2Norm,PreRes)
                     WRITE(88,'(7a)'),'Distance(km)        T_obs(s)         T_ref-iso        T_ref_aa        T_ref &
                             Data_Res           Delta_iso(s)          Res(s) '
                    DO i=1,dall
                        WRITE(88,*) dist(i),obst(i),dsyn(i),RefTaa(i),Tref(i),Tdata(i),fwdTvs(i),resbst(i)
                    ENDDO
                    CLOSE(88)
               ENDIF
               WRITE(*,*)'------------------------||W(Gm-d)||2  finished!--------------------------'
!-----------------------------------------------------------------------------------------------------!
! Statistic for the  data residual.
              mean = sum(resbst(1:dall))/dall
              std_devs = sqrt(sum((resbst(1:dall)-mean)**2)/dall)

                WRITE(66,'(a)')'  '
               WRITE(6,'(i2,a)')  ,iter,'th Inversion Data Fitting.-----After Inversion Data Space------------'
               WRITE(6,*)' Mean Leaving Res (Res_new/Res_org) (%)',PreRes*100
               WRITE(6,'(a,f12.4,a,f12.4,a,f12.4)'),'  (Res_in-T_inv) mean,std_devs and rms of &
                   Res: ',mean*1000,'ms ',1000*std_devs,'ms ',&
                   dnrm2(dall,resbst,1)/sqrt(real(dall))

               WRITE(66,'(i2,a)')  ,iter,'th Inversion Data Fitting.-----After Inversion Data Space------------'
               WRITE(66,*)' Mean Leaving Res (Res_new/Res_org) (%)',PreRes*100
               WRITE(66,'(a,f12.4,a,f12.4,a,f12.4)'),'  (Res_in-T_inv) mean,std_devs and rms of &
                   Res: ',mean*1000,'ms ',1000*std_devs,'ms ',&
                   dnrm2(dall,resbst,1)/sqrt(real(dall))
!-----------------------------------------------------------------------------------------------------!
! Statistic for the inversion space .
               WRITE(*,'(i2,a)')  ,iter,'th Inversion.----------------Model Space-----------------------'
               WRITE(*,'(a,f15.5)') '  RESIDUAL (W(Gm-d)) 2NORM: ',res2Norm
               WRITE(*,'(a,f15.5)') '  RESIDUAL (Gm-d) 2NORM: ',resNSigma2Norm
               WRITE(*,'(a,f15.5)') '  MODEL (Lm) 2NORM: ',Mnorm2
               WRITE(6,'(a,f15.5)') '  MODEL (wLm) 2NORM: ',MwNorm2

               WRITE(66,'(a)')'  '
               WRITE(66,'(i2,a)')  ,iter,'th Inversion.----------------Model Space----------------------'
               WRITE(66,'(a,f15.5)') '  RESIDUAL (W(Gm-d)) 2NORM:    ',res2Norm
               WRITE(66,'(a,f15.5)') '  RESIDUAL (Gm-d) 2NORM: ',resNSigma2Norm
               WRITE(66,'(a,f15.5)') '  MODEL (Lm) 2NORM: ',Mnorm2
               WRITE(66,'(a,f15.5)') '  MODEL  (wLm)  2NORM:         ',MwNorm2

!-----------------------------------------------------------------------------------------------------!
!-----For whole traveltime
!  1. forward calculate the traveltime.
!  2.  forward calculate the phase velocity map.
! Tdata: inversion used data: tobs- tref(iso)
! resbst: residual, without weight(res2Norm with weight)
! Here restT=T_obs-T_forward

 ! WRITE OUT RESIDUAL FOR whole travel time
                WRITE(66,'(a)')'--------------------WHOLE TRAVEL TIME RESIDUAL------------------------------'
                WRITE(id,'(I2.2)') iter
                filename='Traveltime_Result_'//TRIM(id)//'th.dat'
               OPEN(88,file=TRIM(filename))
               resbst=0
               fwdT=0
               IF (iso_inv.eq.0)THEN
                   DO i=1,dall
                        fwdT(i)=dsyn(i)+fwdTaa(i)
                        resbst(i)=obst(i)-fwdT(i)
                   ENDDO
                   WRITE(88,'(7a)'),'Distance(km)       T_obs(s)        T_ref-iso        Res(ref)   &
                           T (Azim)        Res (Azim)         T (aa)'
                   DO i=1,dall
                       WRITE(88,*) dist(i),obst(i),dsyn(i),Tdata(i),fwdT(i),resbst(i),fwdTaa(i)
                   ENDDO
                   FORALL(i=1:dall)
                         RefTaa(i)=fwdTaa(i)
                   END FORALL
              ELSEIF(iso_inv.eq.1)THEN
                  tRcV=0
                  fwdTvs=0
                  GGc=0
                  GGs=0
                  CALL CalSurfGAniso(SenGfile,nx,ny,nz,maxvp,vsf,iw,rw,col,fwdTvs,&
                   GGc,GGs,Lsen_Gsc,dall,rmax,tRcV,&
                   goxd,gozd,dvxd,dvzd,kmaxRc,tRc,periods,depz,minthk,&
                   scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax, nsrc,nrc,nar,Refwritepath)

                  CALL CalAzimTraveltime(nx,ny,nz,maxvp,dall,GGc,GGs,gcf,gsf,fwdTaa)
                  DO i=1,dall
                      fwdT(i)=fwdTvs(i)+fwdTaa(i)
                      resbst(i)=obst(i)-fwdT(i)
                      resbst_iso(i)=obst(i)-fwdTvs(i)
                  ENDDO
                  WRITE(88,'(7a)'),'Distance(km)       T_obs(s)        T_ref-iso        T_ref_aa        T_ref         Res(ref)   &
                           T (Azim)        Res (Azim)       T (iso)       T (aa)'
                  DO i=1,dall
                       WRITE(88,*) dist(i),obst(i),dsyn(i),RefTaa(i),Tref(i),Tdata(i),fwdT(i),resbst(i),fwdTvs(i),fwdTaa(i)
                  ENDDO
                  mean = sum(resbst_iso(1:dall))/dall
                  std_devs = sqrt(sum((resbst_iso(1:dall)-mean)**2)/dall)
                  WRITE(66,'(a,f12.4,a,f12.4,a,f12.4)'),' (Tobs-Tiso) mean,and std_devs  of  whole travel time&
                  residual: ',mean*1000,'ms ',1000*std_devs,'ms '
               ENDIF
               CLOSE(88)
 !-----------------------------------------------------------------------------------------------------!
               mean = sum(resbst(1:dall))/dall
               std_devs = sqrt(sum((resbst(1:dall)-mean)**2)/dall)
               meanAbs=sum(abs(resbst(1:dall)))/dall

               WRITE(66,'(a,f12.4,a,f12.4,a,f12.4)'),' (Tobs-Tiso-Taa) mean,and std_devs  of  whole travel time&
               residual: ',mean*1000,'ms ',1000*std_devs,'ms '


               mean=sum(abs(fwdTaa(1:dall)))/dall
               WRITE(66,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
               WRITE(66,'(a,f12.4,a)')'Azimuthal  Abs Mean Traveltime:',mean,'s'
               WRITE(6,'(a,f12.4)')'Absolute Mean Res(s):',meanAbs
               WRITE(6,'(a,f12.4,a)')'Azimuthal  Abs Mean Traveltime:',mean,'s'
!-----------------------------------------------------------------------------------------------------!
!               deallocate(cbst)
!               WRITE(*,*)'Free cbst'
!               deallocate(dv)
!               WRITE(*,*)'Free dv'
!               deallocate(col)
!                WRITE(*,*)'Free col'
!               deallocate(rw)
!               WRITE(*,*)'Free rw'
!               deallocate(iw)
!               WRITE(*,*)'Free iw'
               WRITE(66,'(i2,a)')  ,iter,'th Inversion finished!'
               WRITE(66,'(a)') ' '
               WRITE(6,'(i2,a)')  ,iter,'th Inversion finished!'
               WRITE(6,'(a)') '  '
           ENDDO ! iteration
!------------------------------------------------------------------------------------------------------!
! output inversion result.
           OPEN(63,file='DSurfTom.inv')
           OPEN(73,file='Gc_Gs_model.inv')
           CALL WriteVsmodel(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,63)
           CALL WriteAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,gcf,gsf,vsf,73)
           CLOSE(63)
           CLOSE(73)
 !-----------------------------------------------------------------------------------------------------!
           OPEN(77,file='period_phaseV_FWD.dat')
           CALL WTPeriodPhaseV(nx,ny,gozd,goxd,dvzd,dvxd,kmaxRc,tRc,tRcV,77)
!-----------------------------------------------------------------------------------------------------!
               WRITE(*,'(a)')'  Begin forward calculate period azimuthal A1, A2.'
               OPEN(42,file='period_Azm_tomo.inv',status='replace',action='WRITE')
               CALL FwdAzimuthalAniMap(nx,ny,nz,maxvp,&
                    goxd,gozd,dvxd,dvzd,kmaxRc,tRc,&
                   gcf,gsf,Lsen_Gsc,tRcV)
!-----------------------------------------------------------------------------------------------------!
! OUTPUT THE VELOCITY MODEL
! NOTE: for the input Vs model is grid model (nz point)
! BUT, the inversion kernel and result is layer (nz-1 layer) model
! for Gc,Gs  the input and result are both layer model, so the depth I set the lower depth as the output
! for Vs, input is grid model, while the inversion result is layered model. For compare withe two Vs models,
!  I set the mid depth between the point as the depth index, and the real model use the average Vs.

            WRITE(*,*),'Program finishes successfully'
            WRITE(66,*),'Program finishes successfully'

            WRITE(*,*),'Output inverted shear velocity model &
                      to Vs_model_Syn.rela and Vs_model_Syn.abs'

            WRITE(66,*),'Output inverted shear velocity model &
                    to Vs_model_Syn.rela and Vs_model_Syn.abs'

        endT=OMP_get_wtime()
        write(*,*)"     All time cost=",endT-startT,"s"
        write(66,*)"     All time cost=",endT-startT,"s"

        CLOSE(10) ! surfaniso.in
        CLOSE(nout) !CLOSE lsmr.txt
        CLOSE(66) !CLOSE surf_tomo.log



        deallocate(obst)
        deallocate(dsyn)
        deallocate(dist)
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
        deallocate(vsf)
        deallocate(sigmaT,Tdata,resbst,fwdT,resSigma,resbst_iso)
        deallocate(fwdTvs,fwdTaa)
        deallocate(RefTaa,Tref)
        deallocate(GVs,GGs,GGc,Lsen_Gsc)

        deallocate(gcf,gsf)
        IF(kmaxRc.gt.0) THEN
        deallocate(tRc)
        ENDIF
        deallocate(tRcV)

        end program

        subroutine WriteVsmodel(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,Idout)
        INTEGER nx,ny,nz
        REAL goxd,gozd
        REAL dvxd,dvzd
        REAL depz(nz)
!        REAL vsAbs(nx-2,ny-2,nz-1), vsRela(nx-2,ny-2,nz-1)
        INTEGER Idout
        REAL:: vsf(nx,ny,nz)
        INTEGER:: k,j,i
        REAL :: vsref
        ! Because the (gozd,goxd) point is NE point of the inversion region.---> NE move 1 point.
        DO k=1,nz
              DO j=1,ny
                    DO i=1,nx
                        WRITE(Idout,'(5f8.4)') gozd+(j-2)*dvzd,goxd-(i-2)*dvxd,depz(k),vsf(i,j,k)
                    ENDDO
                ENDDO
            ENDDO
        end subroutine


        subroutine WriteAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,Gc,Gs,vsf,Idout)
        INTEGER nx,ny,nz
        REAL goxd,gozd
        REAL dvxd,dvzd
        REAL depz(nz)
        REAL Gc(nx-2,ny-2,nz-1),Gs(nx-2,ny-2,nz-1)
        REAL:: vsf(nx,ny,nz)
        INTEGER Idout
        INTEGER:: k,j,i
        real cosTmp,sinTmp
        real vsref
        real*8 :: pi=3.1415926535898
        DO k=1,nz-1
              DO j=1,ny-2
                    DO i=1,nx-2


                         cosTmp=Gc(i,j,k)
                         sinTmp=Gs(i,j,k)
                         AzimAmp=0.5*sqrt(cosTmp**2+sinTmp**2)
                         AzimAng=atan2(sinTmp,cosTmp)/pi*180
                         IF (AzimAng.LT.0.0) AzimAng=AzimAng+360
                         AzimAng=0.5*AzimAng
                         ! format: lon lat  depth(lower) Vs(mid-depth)  Angle  Amp(%/100) Gc/L(%)  Gs/L(%)
                          vsref=(vsf(i+1,j+1,k)+vsf(i+1,j+1,k+1))/2
                        WRITE(Idout,'(8f10.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k+1),vsref,AzimAng,AzimAmp,&
                                   Gc(i,j,k)*100,Gs(i,j,k)*100
                    ENDDO
                ENDDO
            ENDDO
        end subroutine

        subroutine WTPeriodPhaseV(nx,ny,gozd,goxd,dvzd,dvxd,kmaxRc,tRc,pvRc,Idout)
        IMPLICIT NONE
         INTEGER nx,ny
        REAL goxd,gozd
        REAL dvxd,dvzd
        integer kmaxRc
        real*8 tRc(kmaxRc)
        real*8 pvRc((nx-2)*(ny-2),kmaxRc)
        INTEGER Idout
        INTEGER:: tt,k,jj,ii

! NOTE:  pvRc((nx-2)*(ny-2),kmaxRc) ----different pvRc in
             DO tt=1,kmaxRc
                 DO jj=1,ny-2
                     DO ii=1,nx-2
                         WRITE(Idout,'(5f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),pvRc((jj-1)*(nx-2)+ii,tt)
                     ENDDO
                 ENDDO
             ENDDO
             CLOSE(Idout)
        end subroutine

      SUBROUTINE CalRmax(nz,depz,minthk0,rmax)
       IMPLICIT NONE
       INTEGER nz
       REAL depz(nz)
       REAL minthk0,minthk
       INTEGER rmax
       INTEGER k,i
       INTEGER NL
       PARAMETER (NL=200)
       INTEGER nsublay(NL)
       REAL thk
        rmax=0
       DO i=1,nz-1
           thk = depz(i+1)-depz(i)
	   minthk = thk/minthk0
           nsublay(i) = int((thk+1.0e-4)/minthk) + 1
       ENDDO
       rmax=sum(nsublay(1:nz-1))
       END SUBROUTINE
