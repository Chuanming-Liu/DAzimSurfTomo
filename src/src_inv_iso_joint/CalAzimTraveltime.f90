! CalAzimTraveltime:  forward calculate the anisotropic traveltime.
!   used in the Forward program and inversion program.

        SUBROUTINE CalAzimTraveltime(nx,ny,nz,maxvp,dall,GGc,GGs,Gc,Gs,fwdTaa)
        IMPLICIT NONE
! INPUT
        INTEGER,INTENT(IN):: dall,maxvp,nx,ny,nz
        REAL,INTENT(IN):: GGc(dall,maxvp),GGs(dall,maxvp)
        REAL,INTENT(IN):: Gc(nx-2,ny-2,nz-1),Gs(nx-2,ny-2,nz-1)
! OUTPUT
        REAL :: fwdTaa(dall), fwdTgc(dall),fwdTgs(dall)
! PARAMETER
         INTEGER jj,kk,ii,i
         INTEGER nvx,nvz
         REAL GsCol(maxvp),GcCol(maxvp)

        nvx=nx-2
        nvz=ny-2

        DO jj=1,ny-2
        DO kk=1,nx-2
              GcCol( (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=Gc(kk,jj,1:nz-1)
              GsCol( (jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=Gs(kk,jj,1:nz-1)
        ENDDO
        ENDDO


        fwdTgc=MATMUL(GGc,GcCol)
        fwdTgs=MATMUL(GGs,GsCol)
        fwdTaa=0

        DO i=1,dall
                fwdTaa(i)=fwdTgs(i)+fwdTgc(i)
        ENDDO
        END SUBROUTINE

