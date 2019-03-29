! forward period dependent azimuthal anisotropy
! based on the kernel file and input Gc,Gs model
! 1. for true model, with true model and true Gc,Gs;
! 2. for the inversion, with reference model and kernel, and the retrieved Gc,Gs.
!  used in both MainForward.f90 and MainLS.f90
        SUBROUTINE FwdAzimuthalAniMap(nx,ny,nz,maxvp,&
                    goxd,gozd,dvxd,dvzd,kmaxRc,tRc,&
                    gcf,gsf,Lsen_Gsc,tRcV)
        IMPLICIT NONE
! INPUT
	integer nx,ny,nz
        integer maxvp
        real vsf(nx,ny,nz)
        real goxd,gozd,dvxd,dvzd
        integer kmaxRc
        integer writeperiod
        real*8 tRc(kmaxRc)
        real gcf(nx-2,ny-2,nz-1),gsf(nx-2,ny-2,nz-1)
!       character(len=*)SenGfile
       real*8 tRcV((nx-2)*(ny-2),kmaxRc)
        real depz(nz)
        real*4 Lsen_Gsc(nx*ny,kmaxRc,nz-1)
! OUTPUT
! PARAMETER LIST
       real cosRc((nx-2)*(ny-2),kmaxRc),sinRc((nx-2)*(ny-2),kmaxRc)
       INTEGER nvx,nvz
       real cosTmp,sinTmp
       integer ii,jj,kk,j,i,tt,zz,k
       real AzimAmp,AzimAng,AzimAmpRela
        real*8 :: pi=3.1415926535898
        real isoC
!--------------------------------------------------------------------------------!
! LOAD KERNEL
! NOTE: kernel file contains integral kernel.
! AND the kernel file is by ref model which range is nx*ny*(nz-1)-->(nx-2)*(ny-2)
!--------------------------------------------------------------------------------!



! Integral
        nvx=nx-2
        nvz=ny-2
        DO tt=1,kmaxRc
             DO jj=1,ny-2
                  DO ii=1,nx-2
                      cosTmp=0.0
                      sinTmp=0.0
                      DO kk=1,nz-1
                           cosTmp=cosTmp+Lsen_Gsc(jj*(nvx+2)+ii+1,tt,kk)*gcf(ii,jj,kk)
                           sinTmp=sinTmp+Lsen_Gsc(jj*(nvx+2)+ii+1,tt,kk)*gsf(ii,jj,kk)
                      ENDDO
                     cosRc((jj-1)*nvx+ii,tt)=cosTmp
                     sinRc((jj-1)*nvx+ii,tt)=sinTmp
                  ENDDO
             ENDDO
        ENDDO

        !--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
! OUTPUT
! NOTE: cosRc(nvx*nvz,kmaxRc)
       cosTmp=0
       sinTmp=0
!        IF (writeperiod.EQ.1)THEN
             DO tt=1,kmaxRc
                 DO jj=1,ny-2
                     DO ii=1,nx-2
!                         WRITE(41,'(6f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),cosRc((jj-1)*(nx-2)+ii,tt),sinRc((jj-1)*nvx+ii,tt)
                         cosTmp=cosRc((jj-1)*(nx-2)+ii,tt)
                         sinTmp=sinRc((jj-1)*nvx+ii,tt)
                         AzimAmp=sqrt(cosTmp**2+sinTmp**2)
                         isoC=tRcV((jj-1)*(nx-2)+ii,tt)
                         AzimAmpRela=AzimAmp/isoC
                         AzimAng=atan2(sinTmp,cosTmp)/pi*180
                         if (AzimAng.LT.0.0) AzimAng=AzimAng+360
                         AzimAng=0.5*AzimAng
                         !  lon lat period Ciso Angle Amp(%/100) Amp(abs) a1(km/s)  a2(km/s)
                        WRITE(42,'(10f10.5)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),isoC,AzimAng,AzimAmpRela,AzimAmp,&
                        cosTmp,sinTmp
!                       WRITE(42,'(6f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),AzimAmp,AzimAng
                     ENDDO
                 ENDDO
             ENDDO
             close(42)
!        ENDIF

        END SUBROUTINE

