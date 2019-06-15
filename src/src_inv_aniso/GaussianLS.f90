
        subroutine GaussianLS(nz,ny,nx,nar,maxvp,dall ,gozd,goxd,dvzd,dvxd,depz,LcorrXY,LcorrZZ,&
        weightVs,weightGcs,sigmaDeltaVs,sigmaGcs,rw,iw,col,count3,iso_inv)
        implicit none
        integer,intent(in) :: nz,ny,nx,maxvp,dall
        integer :: nar
        real,intent(in) :: goxd,gozd
        real,intent(in) :: dvzd,dvxd
        real,intent(in) :: depz(nz)
        real rw(*)
        integer iw(*),col(*)
        real,intent(in):: LcorrXY,LcorrZZ
        real,intent(in):: weightVs,weightGcs
        INTEGER,INTENT(IN) :: iso_inv
!        real,intent(in):: sigmaDeltaVs(nz-1),sigmaGc(nz-1),sigmaGs(nz-1)
        real,intent(in):: sigmaDeltaVs,sigmaGcs
        integer,intent(out):: count3
!        real,intent(out)guassVs(maxvp,maxvp),guassGc(maxvp,maxvp),guassGs(maxvp,maxvp)

        real :: LcorrXYr
        real :: guassVs(maxvp,maxvp),guassGcs(maxvp,maxvp)
        integer i,j,k,z
        integer nvz,nvx
        integer i1,j1,k1
        integer i2,j2,k2
        real lon1,lon2,lat1,lat2,delta,az,baz

        integer row,column
        real depth1,depth2,DeltaDep
        real LcorrZ1,LcorrZ2

        real exXY,exZZ,exValue
        real twoXYpow2,twoZZpow2
        real covVs,covGcs
        integer INFO
!        real,parameter::ftol=1e-4
        integer count2
        real :: deltakm
        INTEGER::sc
         REAL:: betaVs,betaGcs
        INTEGER:: para1,para2
        REAL :: VsWeight
! gozd,goxd: normal, none colatitude
! delta: Great Circle Arc distance in degrees
!  Latitude of first point (+N, -S) in degrees
! LcorrXY: covert from km to degree
! azdist:
!  maxm =  (nx-2)*(ny-2)*(nz-1)*3 : model parameter number of all model parameters (Vs, Gc and Gs)
!  maxvp = (nx-2)*(ny-2)*(nz-1) : model parameter number of Vs or Gc(Gs)
! Para1: column of model, alphaij ---i main parameter on the cycle, model index on the model
!  Formule: (1/sigma**2)exp[-0.5*(|r-r'|**2)/L**2]
!-------------------------------------------------------!
        nvz=ny-2
        nvx=nx-2
        count3=0
!----------------------------------------------------------------------------------------!
! Step 1. cycle on Vs
IF (iso_inv.eq.1)THEN

        DO para1=1,maxvp
            call CalIndex(para1,nvx,nvz,i1,j1,k1)
             lon1=gozd+(j1-1)*dvzd
             lat1=goxd-(i1-1)*dvxd
             DO para2=1,maxvp
                   IF (para2.NE.para1)THEN
                      call CalIndex(para2,nvx,nvz,i2,j2,k2)
                      lon2=gozd+(j2-1)*dvzd
                      lat2=goxd-(i2-1)*dvxd
                      IF(k1.EQ.k2)THEN   ! SAME Depth(Layer)
                          call azdist(lat1, lon1, lat2, lon2,delta, az, baz)
                          deltakm=delta*111.1949
!                          twoXYpow2=4.*(LcorrXY**2)
                           twoXYpow2=2.*(LcorrXY**2)
                          exXY=(deltakm**2)/twoXYpow2
                          IF (exXY.LT.5)THEN  ! if alpha is not 0
                             betaVs=(1/sigmaDeltaVs)*exp(-exXY)*weightVs
                             !  para Vs
                             count3=count3+1 ! row index
                             rw(nar+1)= -1*betaVs
                             col(nar+1)=para1
                             iw(1+nar+1)=dall+count3
                             ! para2
                             rw(nar+2)= betaVs
                             col(nar+2)=para2
                             iw(1+nar+2)=dall+count3
                             nar=nar+2  ! G index
                          ENDIF
                      ELSE IF((i1.EQ.i2).AND.(j1.EQ.j2) )THEN  ! SAME GRID 2D, but not same depth(not same point)
                          depth1=(depz(k1)+depz(k1+1))/2
                          depth2=(depz(k2)+depz(k2+1))/2
                          DeltaDep=abs(depth1-depth2)
                          twoZZpow2=2.*(LcorrZZ**2)
                          exZZ=(DeltaDep**2)/twoZZpow2

                          IF (exZZ.LT.5)THEN  ! if alpha is not 0
                             betaVs=(1/sigmaDeltaVs)*exp(-exZZ)*weightVs
                             !  para1
                             count3=count3+1 ! row index
                             rw(nar+1)=-1*betaVs
                             col(nar+1)=para1
                             iw(1+nar+1)=dall+count3
                             ! para2
                             rw(nar+2)= betaVs
                             col(nar+2)=para2
                             iw(1+nar+2)=dall+count3
                             nar=nar+2  ! G number index
                          ENDIF
                      ENDIF ! (k1.EQ.k2)
                   ENDIF !para1=para2
             ENDDO ! para2
        ENDDO !para1
!-------------------------------------------------------------------------------------------------!
! Gcs-- First cycle is Gc, second cycle is Gs
ELSEIF (iso_inv.eq.0)THEN
        DO sc=1,2
            DO para1=1,maxvp
               call CalIndex(para1,nvx,nvz,i1,j1,k1)
               lon1=gozd+(j1-1)*dvzd
               lat1=goxd-(i1-1)*dvxd
               DO para2=1,maxvp
                   IF (para2.NE.para1)THEN
                        call CalIndex(para2,nvx,nvz,i2,j2,k2)
                         lon2=gozd+(j2-1)*dvzd
                         lat2=goxd-(i2-1)*dvxd
                         IF(k1.EQ.k2)THEN   ! SAME Depth(Layer)
                             call azdist(lat1, lon1, lat2, lon2,delta, az, baz)
                             deltakm=delta*111.1949
!                             twoXYpow2=4.*(LcorrXY**2)
                               twoXYpow2=2.*(LcorrXY**2)
                             exXY=(deltakm**2)/twoXYpow2
                             IF (exXY.LT.5)THEN  ! if alpha is not 0
                                betaGcs=(1/sigmaGcs)*exp(-exXY)*weightGcs
                                !  para
                                count3=count3+1 ! row index
                                rw(nar+1)=betaGcs
                                col(nar+1)=para1+(sc-1)*maxvp
                                iw(1+nar+1)=dall+count3
                             ! para2
                                rw(nar+2)= -1*betaGcs
                                col(nar+2)=para2+(sc-1)*maxvp
                                iw(1+nar+2)=dall+count3
                                nar=nar+2  ! G index
                          ENDIF
                      ELSE IF((i1.EQ.i2).AND.(j1.EQ.j2)) THEN! SAME GRID 2D, but not same depth(not same point)
                          depth1=(depz(k1)+depz(k1+1))/2
                          depth2=(depz(k2)+depz(k2+1))/2
                          DeltaDep=abs(depth1-depth2)
                          twoZZpow2=2.*(LcorrZZ**2)
                          exZZ=(DeltaDep**2)/twoZZpow2
                          IF (exZZ.LT.5)THEN  ! if alpha is not 0
                             betaGcs=(1/sigmaGcs)*exp(-exZZ)*weightGcs
                             !  para
                             count3=count3+1 ! row index
                             rw(nar+1)=betaGcs
                             col(nar+1)=para1+(sc-1)*maxvp
                             iw(1+nar+1)=dall+count3
                             ! para2
                             rw(nar+2)= -1*betaGcs
                             col(nar+2)=para2+(sc-1)*maxvp
                             iw(1+nar+2)=dall+count3
                             nar=nar+2  ! G index
                          ENDIF
                      ENDIF ! (k1.EQ.k2)
                   ENDIF !para1=para2
             ENDDO ! para2
        ENDDO !para1
      ENDDO ! sc=1,2rese
ENDIF

       end subroutine GaussianLS

       subroutine CalIndex(para,nvx,nvz,i,j,k)
       implicit none
       integer,intent(in):: para
       integer,intent(in):: nvx,nvz
       integer,intent(out):: i,j,k
       integer:: res

       k=int(para/(nvz*nvx))
       k=k+1
       res=para-(k-1)*(nvz*nvx)
       if (res.eq.0)then
            k=k-1
            res=para-(k-1)*(nvz*nvx)
       endif
       j=int(res/(nvx))
       j=j+1
       i=res-(j-1)*nvx
        if (i.eq.0)then
            j=j-1
            i=res-(j-1)*nvx
         endif
         end subroutine CalIndex


       subroutine corlengthDepth(depth,corlength)
       implicit none
       real:: depth,corlength
       real:: depthMin=0.0,depthMax=100
       real,parameter::ftol=1e-4
       real:: corlengthMin=30,corlengthMax=150
!---------------------------------------------
       if (abs(depth-depthMin).le.ftol ) corlength=corlengthMin
       if (depth >= depthMax ) corlength=corlengthMax
       if ( (depth> depthMin) .and. (depth<depthMax) ) then
            corlength=corlengthMin+(depth-depthMin)*(corlengthMax-corlengthMin)/(depthMax-depthMin)
       end if
       corlength=100
       end subroutine corlengthDepth

