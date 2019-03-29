! CODE FOR SURFACE WAVE  ANISOTROPY TOMOGRAPHY USING DISPERSION MEASUREMENT
! MAKE ANISOTROPIC KERNEL dc/dL  dc/dA
! VERSION: 1.0
! AUTHOR:
! CHUANMING LIU. chuanmingliu@foxmail.com
! HISTORY:
!        2016/10/22

       PROGRAM Read_DKernel
       IMPLICIT NONE
! VARIABLE DEFINE
            character dummy*40
            character*500 refmodel_file

            integer nx,ny,nz
            integer kmax,kmaxRc
            integer checkstat
            integer:: i,j,k,z,ff
            real*8,dimension(:,:,:),allocatable::LonSen_dcdL,LonSen_dcdA
            real*8,dimension(:,:,:),allocatable::LonSen_dcdC,LonSen_dcdF
            real*8,dimension(:,:,:),allocatable::sen_dcdL,sen_dcdA
            real*8,dimension(:,:,:),allocatable::sen_dcdC,sen_dcdF
            integer:: ny_clc,lay,rmax
            character Lonlist*100
            character Lonfile*100
! Main Variables List
! depz(nz): model depth array, contains 0 km
! vsf(nx,ny,nz): input model, contains the out cycle grid which will participate in inversion, but need to cal sensitivity kernel.
!-----------------------------------------
! OUTPUT PROGRAM INFOMATION
             write(*,*)
             write(*,*),'                  Read DKernel'
             write(*,*),'PLEASE contact Chuanming Liu '
             write(*,*)

! 1.0 READ INPUT FILE
            write(*,*) 'model dimension:nx,ny,nz:'
            read(*,*) nx,ny,nz
            write(*,'(3i5)') nx,ny,nz
            write(*,*)'Rayleigh data period number:'
            read(*,*) kmaxRc
            write(*,*)kmaxRc
            write(*,*)'Longitude kernel list file:'
            read(*,*) Lonlist
            write(*,*)Lonlist
            open(34,file=Lonlist,status='old',action='read')

            do ff=1,ny
               read(34,*)Lonfile
               open(35,file=Lonfile,status='old',action='read')
               read(35,*)ny_clc
               read(35,*)rmax

               IF (ff.eq.1) THEN
                 allocate(LonSen_dcdL(nx,kmaxRc,rmax),LonSen_dcdA(nx,kmaxRc,rmax), stat=checkstat)
                 allocate(LonSen_dcdC(nx,kmaxRc,rmax),LonSen_dcdF(nx,kmaxRc,rmax), stat=checkstat)
                 allocate(sen_dcdL(ny*nx,kmaxRc,rmax),sen_dcdA(ny*nx,kmaxRc,rmax), stat=checkstat)
                 allocate(sen_dcdC(ny*nx,kmaxRc,rmax),sen_dcdF(ny*nx,kmaxRc,rmax), stat=checkstat)
               ENDIF

               read(35,*)(((LonSen_dcdL(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
               read(35,*)(((LonSen_dcdA(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
               read(35,*)(((LonSen_dcdC(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
               read(35,*)(((LonSen_dcdF(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
               close(unit=35)

               sen_dcdL((ny_clc-1)*nx+1:(ny_clc)*nx,1:kmaxRc,1:rmax)=LonSen_dcdL(1:nx,1:kmaxRc,1:rmax)
               sen_dcdA((ny_clc-1)*nx+1:(ny_clc)*nx,1:kmaxRc,1:rmax)=LonSen_dcdA(1:nx,1:kmaxRc,1:rmax)
               sen_dcdC((ny_clc-1)*nx+1:(ny_clc)*nx,1:kmaxRc,1:rmax)=LonSen_dcdC(1:nx,1:kmaxRc,1:rmax)
               sen_dcdF((ny_clc-1)*nx+1:(ny_clc)*nx,1:kmaxRc,1:rmax)=LonSen_dcdF(1:nx,1:kmaxRc,1:rmax)
            enddo
            close(unit=34)

            open(37,file='Sen_dcdL_dA_dC_dF.dat',action='write')
            DO j=1,nx*ny
               write(37,*)((sen_dcdL(j,i,z),z=1,rmax),i=1,kmaxRc)
            ENDDO
            DO j=1,nx*ny
               write(37,*)((sen_dcdA(j,i,z),z=1,rmax),i=1,kmaxRc)
            ENDDO
             DO j=1,nx*ny
               write(37,*)((sen_dcdC(j,i,z),z=1,rmax),i=1,kmaxRc)
            ENDDO
            DO j=1,nx*ny
               write(37,*)((sen_dcdF(j,i,z),z=1,rmax),i=1,kmaxRc)
            ENDDO
            close(unit=37)

            !close surf_tomo.log
            write(*,*)'DK reader is over.^-^'
            deallocate(LonSen_dcdL,LonSen_dcdA)
            deallocate(LonSen_dcdC,LonSen_dcdF)
            deallocate(sen_dcdL,sen_dcdA)
            deallocate(sen_dcdC,sen_dcdF)


       END PROGRAM Read_DKernel

