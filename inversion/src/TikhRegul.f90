! ADDING REGULARIZATION TERM



        SUBROUTINE TikhonovRegularization(nx,ny,nz,maxvp,dall,nar,rw,iw,col,count3,iso_inv,weightGcs,weightVs)
        IMPLICIT NONE
        INTEGER nx,ny,nz,maxvp
        INTEGER dall,nar
        REAL rw(*)
        INTEGER iw(*),col(*)
        INTEGER iso_inv,count3
        REAL weightGcs, weightVs
        REAL weight

        INTEGER i,j,k
        INTEGER nvx,nvz,sc
   !-------------------------------------------------------!
        nvz=ny-2
        nvx=nx-2
        count3=0
!----------------------------------------------------------!
        IF (iso_inv.eq.1)THEN
                weight=weightVs
                do k=1,nz-1
                 do j=1,nvz
                 do i=1,nvx
                 if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
                 count3=count3+1
                 col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
                 rw(nar+1)=2.0*weight
                 iw(1+nar+1)=dall+count3

                 nar=nar+1
                 else
                 count3=count3+1
                 col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
                 rw(nar+1)=6.0*weight
                 iw(1+nar+1)=dall+count3
                 rw(nar+2)=-1.0*weight
                 iw(1+nar+2)=dall+count3
                 col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1
                 rw(nar+3)=-1.0*weight
                 iw(1+nar+3)=dall+count3
                 col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1
                 rw(nar+4)=-1.0*weight
                 iw(1+nar+4)=dall+count3
                 col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i
                 rw(nar+5)=-1.0*weight
                 iw(1+nar+5)=dall+count3
                 col(nar+5)=(k-1)*nvz*nvx+j*nvx+i
                 rw(nar+6)=-1.0*weight
                 iw(1+nar+6)=dall+count3
                 col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i
                 rw(nar+7)=-1.0*weight
                 iw(1+nar+7)=dall+count3
                 col(nar+7)=k*nvz*nvx+(j-1)*nvx+i

                 nar=nar+7
                 endif
                 enddo
                 enddo
                 enddo
        ELSEIF (iso_inv.eq.0)THEN
              weight=weightGcs
              DO sc=1,2
                 do k=1,nz-1
                 do j=1,nvz
                 do i=1,nvx
                 if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
                 count3=count3+1
                 col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
                 rw(nar+1)=2.0*weight
                 iw(1+nar+1)=dall+count3

                 nar=nar+1
                 else
                 count3=count3+1
                 col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
                 rw(nar+1)=6.0*weight
                 iw(1+nar+1)=dall+count3
                 rw(nar+2)=-1.0*weight
                 iw(1+nar+2)=dall+count3
                 col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1+(sc-1)*maxvp
                 rw(nar+3)=-1.0*weight
                 iw(1+nar+3)=dall+count3
                 col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1+(sc-1)*maxvp
                 rw(nar+4)=-1.0*weight
                 iw(1+nar+4)=dall+count3
                 col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i+(sc-1)*maxvp
                 rw(nar+5)=-1.0*weight
                 iw(1+nar+5)=dall+count3
                 col(nar+5)=(k-1)*nvz*nvx+j*nvx+i+(sc-1)*maxvp
                 rw(nar+6)=-1.0*weight
                 iw(1+nar+6)=dall+count3
                 col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
                 rw(nar+7)=-1.0*weight
                 iw(1+nar+7)=dall+count3
                 col(nar+7)=k*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp

                 nar=nar+7
                 endif
                 enddo
                 enddo
                 enddo
              ENDDO
           ENDIF

        END SUBROUTINE






