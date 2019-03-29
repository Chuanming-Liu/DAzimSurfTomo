! CODE FOR SURFACE WAVE  ANISOTROPY TOMOGRAPHY USING DISPERSION MEASUREMENT
! MAKE ANISOTROPIC KERNEL dc/dL  dc/dA
! VERSION: 1.0
! AUTHOR:
! CHUANMING LIU. chuanmingliu@foxmail.com
! HISTORY:
!        2016/10/22

       PROGRAM Aniso_kernel
       IMPLICIT NONE
! VARIABLE DEFINE
           character inputfile*80
            character mldfile*400
            character logfile*100
            logical ex
            character dummy*40
            character*500 refmodel_file

            integer nx,ny,nz
            real goxd,gozd
            real dvxd,dvzd
            real minthk
            integer kmax,kmaxRc,kmaxRg,kmaxLc,kmaxLg
            real*8,dimension(:),allocatable:: tRc,tRg,tLc,tLg
            real,dimension(:),allocatable:: depz
            real,parameter :: pi=3.1415926535898
            integer checkstat
            integer:: i,j,k
            integer,dimension(:,:),allocatable::periods
            real,dimension(:,:,:),allocatable::vsf
            real*8,dimension(:,:,:),allocatable::sen_LRc,sen_ARc
            integer:: ny_clc,rmax
! Main Variables List
! depz(nz): model depth array, contains 0 km
! vsf(nx,ny,nz): input model, contains the out cycle grid which will participate in inversion, but need to cal sensitivity kernel.
!-----------------------------------------
! OUTPUT PROGRAM INFOMATION
             write(*,*)
             write(*,*),'                         ANISO KERNEL'
             write(*,*),'PLEASE contact Chuanming Liu &
                 (chuanmingliu@foxmail.com) if you find any bug.'
             write(*,*)

! 1.0 READ INPUT FILE
            write(*,*) 'input file [AniKernel.in(default)]:'
            write(*,*)'initial model:'
            read(*,'(A)')mldfile
            write(*,*)mldfile
            write(*,*) 'model dimension:nx,ny,nz:'
            read(*,*) nx,ny,nz
            write(*,'(3i5)') nx,ny,nz
            write(*,*)'Longitude cycle num:'
            read(*,*) ny_clc
            write(*,*)ny_clc
            write(*,*)'minthk= layerthick / sublayernum, can be 1, 2 ,3'
            read(*,*)minthk
            write(*,*)minthk
            write(*,*)'Rayleigh data period number:'
            read(*,*) kmaxRc
            write(*,*)kmaxRc

            inputfile='AniKernel'
            write(logfile,'(a,a)')trim(inputfile),'.log'
            open(66,file=logfile)
            write(66,*)
            write(66,*),'                      ANISO KERNEL'
            write(66,*),'PLEASE contact   Chuanming Liu &
                  (chuanmingliu@foxmail.com) if you find any bug'
            write(66,*)'initial model:MOD'
            write(66,'(a80)')mldfile
            write(66,*) 'model dimension:nx,ny,nz'
            write(66,'(3i5)') nx,ny,nz
            write(66,*) 'model dimension longitude ny'
            write(66,*)ny_clc

            if(kmaxRc.gt.0)then
               allocate(tRc(kmaxRc),&
                stat=checkstat)
               if (checkstat > 0) stop 'error allocating RP'
            read(*,*)(tRc(i),i=1,kmaxRc)
            write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
            write(*,'(50f6.2)')(tRc(i),i=1,kmaxRc)
            write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
            write(66,'(50f6.2)')(tRc(i),i=1,kmaxRc)
            endif

            write(*,*)'Mineos reference model:'
            read(*,'(A)')refmodel_file
            write(*,'(a80)')refmodel_file
            write(66,*)'Mineos reference model:'
            write(66,'(a80)')refmodel_file

! 2.0 MEASUREMENTS STATISTICS AND READ INITIAL MODEL
            allocate(depz(nz), stat=checkstat)
            allocate(vsf(nx,ny,nz), stat=checkstat)
            open(10,file=mldfile,status='old')
            read(10,*) (depz(i),i=1,nz)
            do k = 1,nz
            do j = 1,ny
            read(10,*)(vsf(i,j,k),i=1,nx)
            enddo
            enddo
            close(10)
            write(*,*) 'grid points in depth direction:(km)'
            write(*,'(50f6.2)') depz

           Call CalRmax(nz,depz,minthk,rmax)
! 3.0  CALCULATE SENSITIVITY KERNEL BASED ON THE INPUT MODEL
            CALL  CalSenKernel(nx,ny,nz,rmax, vsf,  kmaxRc, tRc,&
             depz,minthk,refmodel_file,ny_clc)
            !close surf_tomo.log
            write(*,*)'DK is over.^-^'
            write(66,*)'DK is over.'
            close(66)
       END PROGRAM Aniso_kernel


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


