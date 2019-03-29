!--------------------------------------------------------------------------------------------!
!  target:
! 1.  interpolate the model to grid models
! 2.  add the perm model to the core
! 3.  output Mineos format model and reference depth files.
!--------------------------------------------------------------------------------------------!
      subroutine  mineos_model_output(refmodel_file,outmdl,&
       rmax,rdep,rvp,rvs,rrho,rthk,kmaxRc,tRc)
       implicit none

      integer,parameter:: NLAY=200,NL=200
      integer,parameter:: NP=512
!     length of the grid model
      integer,parameter:: N_GRID=400
!     length of the interpolate model
      integer,parameter:: NINTP=50000
      real:: epslon
      real*8 ::  Radius=6371.00
      real*4 :: mks=1000
!     input variable
      character(len=*) refmodel_file
      character(len=*) outmdl
      integer:: rmax
      real :: rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
      integer kmaxRc
      real*8 tRc(kmaxRc)
!    SETTING
      real:: thresholdT,minT, limT
      real :: thDep
      real*4 :: dz(2)
      real*4 :: depth,move
      real*4 :: thick,minthk,newthk
      integer:: kk,nsublay
      real*4:: x,x1,x2,y
      integer:: i,j,p,z,k
!     variable in temp model
      real*4 model(NLAY,4)  ! original model
      real*4 model_comb(NLAY,4) ! combine the same layer
      real*4 model_grid(N_GRID,4) ! covert the layer model to grid model
      real*4 model_int(NINTP,4)  ! the interpolated model filtered with depth interval dz km
      real*4  model_up(NINTP,4)  ! the m-kg-s model  in depth descending order,radius increasing.
      real*4 model_ref(NINTP,4)  ! the reference perm model , 2km interval
      real*4  temp_model(4)
      real*4  mod_margin(2,4) ! the margin model
      real*4  model_fl(NINTP,4) ! the output model  model_ref+model_mid+model_up
      real*4  model_mid(NLAY,4) ! the middle model between model_ref and model_up

      real*4  M_mdl(NINTP,9)
!     layer or grid
      integer:: maxlayer,Nlayer,DNlayer,Ngrid
      REAL*4 qkappa(NINTP),qshear(NINTP)
      REAL*4 Qscan(NINTP,4)
      REAL*4 Qp(NINTP),Qs(NINTP)
      real*4 depth_out(NLAY)
!     real(kind=8) model_ref(ndata,4)
      integer n_ref,nic_ref,noc_ref
      integer mid_length
      real*4 :: up_margin,dw_margin
      integer:: Ntrun,locate,Nmid,Nknot

      real*4 :: eta
      real*8 tref
      integer:: ifanis,ifdeck
      integer:: nic,noc,lrec
      integer:: pt,indexq,moho_ind
      LOGICAL:: isexist

!----------------------------------------------------------------------
! VARIABLE LIST
! INPUT LIST
!  rmax: number of layers in the fined layered model
!  rvp, rvs, rrho,: the refined layered velocity model
!  Note: dep contains 0km, but rdep does not contain 0Km
!  rdep:depth of lower layer in refined layered model
!  rthk: thickness of lower layer in refined layered model
!  SETING
!  thresholdT: threshold period (s)
!  limT: limited period (s)
!  minT: minimum value of tRc (s)
!  limT----thresholdT
!  thDep: interpolated boundary (km)
!  SW: switch of whether to interpolate in upper layer(<thDep), 1: need, 2: not need
!  dz(1): <=thDep, shollow interpolation level(km); dz(2)>theDep, deeper interpolation level
!   move: the staggered distance between the interface, when covert layer model to grid model.
!  TEMPE LIST
!  maxlayer: number of layer of model, last layer is NOT half-space(deleteed half-space)
!  model: initial model, [thickness,vp,vs,rho]
!  model_comb: [thick,vp,vs,rho]--- layer model
!  model_grid:[depth,vp,vs,rho]--- grid model
!  Nlayer: number of layer of model_comb
!  DNlayer=Nlayer*2: number of grid of model_grid
!  Grid number:
! Nknot=Ntrun+Nmid+Ngrid
!  Ntrun: for the lower part of the model use the reference model directly
!  Nmid: length of the transition zone.
!  Ngrid: length of the interpolated grid model.
!  Nup=Nmid+Ngrid !\ upper part of model
!
!   mid_length: length of mid layer. Point Number.
!       outmdl='Mieos_process_model.card'
!       out_depth='Depth_segment.dat'
!---------------------------------------------------------------------
! SETTING
       thresholdT=2
       limT=1
       minT=minval(tRc(1:kmaxRc))
       thDep=5
        dz=0
       IF (minT<limT)THEN
         write(*,*)'Error! minT must >',thresholdT
         stop
       ELSE IF (minT>=thresholdT) THEN
               dz(1)=0.5
               dz(2)=0.5
       ELSE IF ((minT<thresholdT).and.(minT>=limT))THEN
                dz(1)=0.1
                dz(2)=0.5
       ENDIF
!---------------------------------------------------------------------
!   step1 . combine the same layer
!   model:[thick,vp,vs,rho]---layer model, Nlayer
!   For refined layered model, should not contain same layer
!   which means maxlayer should equal to Nlayer
!   and in fact do not need to combine the same layer.
!---------------------------------------------------------------------
       epslon = 1.0d-6
       depth_out(1)=0.0
       maxlayer=rmax
       Do i=1, maxlayer
               model(i,1)=rthk(i)
               model(i,2)=rvp(i)
               model(i,3)=rvs(i)
               model(i,4)=rrho(i)
               depth_out(i+1)=depth_out(i)+rthk(i)
       Enddo

!  model_comb=[thick,vp,vs,rho]--- layer model
       model_comb(1,1:4)=model(1,1:4)
       Nlayer=1
       Do i=2, maxlayer
             temp_model(1:4)=model(i,1:4)
             if (temp_model(2).eq.model_comb(Nlayer,2)) then
                model_comb(Nlayer,1)=model_comb(Nlayer,1)+temp_model(1)
             else
                Nlayer=Nlayer+1
                model_comb(Nlayer,1:4)=temp_model(1:4)
                endif
       Enddo
!-----------------------------------------------------------------------------
!  step 2. interpolate the interface, change the layer model to grid model
!  model_grid=[depth,vp,vs,rho] grid model, DNlayer:number of model_grid knot
!  NOTE: 1. move: mainly limit the too thick layer in depth;
!--------------------------------------------------------------------------------
       model_grid(1,1)=0
       model_grid(1,2:4)=model_comb(1,2:4)
      ! cycle on the interface, depth means the last interface depth
       depth=0;
       Do i=1,Nlayer-1   !  i= index of boundary.
         move=model_comb(i,1)/10
         if (depth.LE.thDep.and.  move.GT.dz(1))move=dz(1)
         if (depth.GT.thDep.and.move.GT.dz(2))move=dz(2)

         depth=depth+model_comb(i,1) ! depth of the i th interface
         model_grid(i*2,1)=depth-move/2       ! upper interface
         model_grid(i*2,2:4)=model_comb(i,2:4)
         model_grid(i*2+1,1)=depth+move/2   ! lower interface
         model_grid(i*2+1,2:4)=model_comb(i+1,2:4)
       ENDDO
       DNlayer=2*Nlayer
       model_grid(DNlayer,1)=depth+model_comb(Nlayer,1)
       model_grid(DNlayer,2:4)=model_comb(Nlayer,2:4)

!       OPEN(unit=55,FIlE='model_grid.dat',ACTION='WRITE')
!       WRITE(unit=55,FMT='(4F8.3)')((model_grid(i,j),j=1,4),i=1,DNlayer)
!       CLOSE(UNIT=55)

!-----------------------------------------------------------------------------
!  step 3.  TEST INTERPOLATION
!  IF the thickness between two grid is great than dz
!  interpolate between the two grids.
!  model_int=[depth,vp,vs,rho]-----> interpolated model
!  grid num=
!-----------------------------------------------------------------------------
       model_int(1,1:4)=model_grid(1,1:4) !depth=0
       kk=1
       x=0
       do i=2,DNlayer
           x1=model_grid(i-1,1)  ! upper grid, record
           x2=model_grid(i,1)      ! lower grid
           thick=x2-x1
           if (x1.LE.thDep)then
               minthk=dz(1)
           else
                minthk=dz(2)
           endif
           nsublay=nint((thick+1.0e-4)/minthk)
           if (nsublay.EQ.0) nsublay=1
           newthk=thick/nsublay
           do j=1,nsublay
              kk=kk+1
              x=x+newthk
              model_int(kk,1)=x
              do p=2,4
                 call linear_intp_dep(x1,model_grid(i-1,p),x2,model_grid(i,p),x,y)
                  model_int(kk,p)=y
              enddo
           enddo
        enddo
        Ngrid=kk

!       OPEN(unit=56,FIlE='model_int.dat',ACTION='WRITE')
!       WRITE(unit=56,FMT='(4F8.3)')((model_int(i,j),j=1,4),i=1,Ngrid)
!       CLOSE(UNIT=56)

!        convert kgs to mks
!        For model_grid k-g-s unit (km-GM/CC(g/cm^3)-s)
!        For mineos model format : m-kg-s (m-kg/m^3-s)

!        model_up=[radius,vp,vs,rho]: upper model, depth descending order.
        Do z=1,Ngrid
            model_up(z,1)=(Radius-model_int(Ngrid-z+1,1))*mks
            model_up(z,2:4)=model_int(Ngrid-z+1,2:4)*mks
        Enddo

!       OPEN(unit=55,FIlE='input_model.dat',ACTION='WRITE')
!       WRITE(unit=55,FMT='(4F12.2)')((model_up(i,j),j=1,4),i=1,Ngrid)
!       CLOSE(UNIT=55)

!        write(*,'(F9.2)')(model_up(i,1),i=1,Ngrid)
!-----------------------------------------------------------------------------------------------
!        Step 4. load the Mineos card for the mental and core model
!        model_ref=[radius,vp,vs,rho,Qp,Qs]
       INQUIRE(FILE=refmodel_file,EXIST=isexist)
       IF ( isexist ) then
           CALL read_reference_model(refmodel_file,&
            model_ref, n_ref,nic_ref,noc_ref,qkappa,qshear)
       Else
            write(*,*)refmodel_file
            STOP 'SUB:mineos_model_output :Can not find the reference model file!'
       Endif
!        Step 5. combine the ref_perm model (lower part core mantle) and upper model
!       mid_length=10---10*dz(2): 5km gap depth mid
       mid_length=100 ! 50km gap.
       up_margin=model_up(1,1)-mid_length*dz(2)*mks
       dw_margin=model_up(1,1)
       mod_margin(2,1:4)=model_up(1,1:4) !down margin
       Ntrun=locate(model_ref(1:n_ref,1),n_ref,up_margin)
       mod_margin(1,1:4)=model_ref(Ntrun,1:4)  ! up
       up_margin=model_ref(Ntrun,1)
       Nmid=ceiling((dw_margin-up_margin)/(dz(2)*mks))-1

!     model_mid=[radius,vp,vs,rho]-> mid trans zone, depth descending order.
       Do i=1,Nmid
          model_mid(i,1)=up_margin+dz(2)*mks*i
       Enddo
       Do k=2,4
          Do i=1,Nmid
              call linear_intp_dep(mod_margin(1,1),mod_margin(1,k),mod_margin(2,1),&
              mod_margin(2,k),model_mid(i,1),model_mid(i,k))
          Enddo
       Enddo
!    index of model
        Nknot=Ntrun+Nmid+Ngrid !Nknot is the number of model knots
!       Ntrun: for the lower part of the model use the reference model directly
!       Nmid: length of the transition zone.
!       Ngrid: length of the interpolated grid model.
!       Nup=Nmid+Ngrid !\ upper part of model

!      model_fl=[radius,vp,vs,rho]--> complete grid model, depth descending order.
       model_fl(1:Ntrun,1:4)=model_ref(1:Ntrun,1:4)
       model_fl(Ntrun+1:Ntrun+Nmid,1:4)=model_mid(1:Nmid,1:4)
       model_fl(Ntrun+Nmid+1:Ntrun+Nmid+Ngrid,1:4)= model_up(1:Ngrid,1:4)
!%----------------------------------------------------------------------------------------------------------------------------------%!
!     Q value part
!    Q scan fault: For the Q value only change when counter discontinuity( double radius)
!    Qscan =[RadiusS,RadiusE,Qp,Qs]
!%----------------------------------------------------------------------------------------------------------------------------------%!
       Qp(1:Ntrun)=qkappa(1:Ntrun) !\ ref model
       Qs(1:Ntrun)=qshear(1:Ntrun)


       Qscan(1,1)=model_ref(1,1)
       Qscan(1,3)=qkappa(1)
       Qscan(1,4)=qshear(1)
        k=1
       DO i=2,n_ref
           IF (model_ref(i,1).EQ.model_ref(i-1,1) )THEN
                Qscan(k,2)=model_ref(i,1)
                k=k+1
                Qscan(k,1)=model_ref(i,1)
                Qscan(k,3)=qkappa(i)
                Qscan(k,4)=qshear(i)
             ENDIF
        ENDDO
        Qscan(k,2)=model_ref(n_ref,1)

        DO i=Ntrun+1,Nknot
           DO j=1,k
               IF (Qscan(j,2) >= model_fl(i,1)) THEN
                   Qp(i)=Qscan(j,3)
                   Qs(i)=Qscan(j,4)
                   EXIT
              ENDIF
            ENDDO
        ENDDO

!%----------------------------------------------------------------------------------------------------------------------------------%!

!      Step 6. output model for Mineos and depth  file
!      model_fl=[radius,vp,vs,rho,Qp,Qs]
!      M_mdl: [radius, rho, vpv, vsv,qkapp,qshear,vph,vsh,eta
!
       eta=1.0
       M_mdl(1:Nknot,1)=model_fl(1:Nknot,1)
       M_mdl(1:Nknot,2)=model_fl(1:Nknot,4)
       M_mdl(1:Nknot,3:4)=model_fl(1:Nknot,2:3) ! vpv,vsv
       M_mdl(1:Nknot,5)=Qp(1:Nknot)
       M_mdl(1:Nknot,6)=Qs(1:Nknot)
       M_mdl(1:Nknot,7:8)=model_fl(1:Nknot,2:3) ! vph,vsh
       M_mdl(1:Nknot,9)=eta
       ifanis=1 !   ifanis=1 for anisotropic model, 0 for isotropic
!       tref= -1 !tref=ref period(secs) of model for dispersion correction.. If tref ≤ 0, no correction is made
       call calculate_tref(kmaxRc,tRc,tref)
       ifdeck= 1 !   ifdeck=1 for card deck model, 0 for polynomial model.
       nic=nic_ref ! nic is the index of the solid side of the inner core  boundary (ICB).
       noc=noc_ref  !noc is the index of the fluid side of the mantle core boundary (MCB).
       lrec=0  ! Irec=1,0---saveall=0
       OPEN(UNIT=41,FILE=outmdl,ACTION='WRITE')
       write(UNIT=41,FMT='(A)')'Isotropic model'
       write(UNIT=41,FMT='(I3,f6.2,I3)')ifanis,tref,ifdeck
       write(UNIT=41,FMT='(4I8)')Nknot,nic,noc,lrec
       write(UNIT=41,FMT='(f9.1,3f10.2,2f10.1,2f10.2,f10.2)')((M_mdl(i,j),j=1,9),i=1,Nknot)
       Close(UNIT=41)

      End Subroutine mineos_model_output



!----------------------------------------------------------------------------------------!
! subroutine read_reference_model load the model file of mineos format
!-----------------------------------------------------------------------------------------!
    subroutine read_reference_model(refmodel_file,model_ref,&
        n_ref,nic_ref,noc_ref,qkappa,qshear)
!        include'mineos.inc'
        implicit none
	character(len=500),intent(in):: refmodel_file
	character(len=40) ititle
	integer ifanis,tref,ifdeck
        integer,intent(out) :: n_ref,nic_ref,noc_ref
        Integer Irec
	integer,parameter:: ndata=50000
	real(kind=8) r(ndata),rho(ndata),vpv(ndata),vsv(ndata),&
	           vph(ndata),vsh(ndata),eta(ndata)
        REAL(kind=4),INTENT(OUT):: qkappa(ndata),qshear(ndata)
           integer ::iin=41
           real(kind=4),intent(out):: model_ref(ndata,4)
           integer:: n,i

       open(iin,file=refmodel_file,status='old',form='formatted')
	   read(iin,*) ititle
       read(iin,*) ifanis,tref,ifdeck
       read(iin,*) n_ref,nic_ref,noc_ref,Irec
!        *** card deck model ***
       read(iin,*) (r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i),&
               vph(i),vsh(i),eta(i),i=1,n_ref)
! 1055     format(f8.0,3f9.2,2f9.1,2f9.2,f9.5)
!	     read(iin,1055) (r(i),rho(i),vpv(i),vsv(i),
!     1    qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
!    1055	format(f8.0,3f9.2,2f9.1,2f9.2,f9.5)
       close(iin)
       Do i=1,n_ref
          model_ref(i,1)=r(i)  ! radius
          model_ref(i,2)=vpv(i)  ! vp
          model_ref(i,3)=vsv(i)  ! vs
          model_ref(i,4)=rho(i)  ! rho
       Enddo
      end subroutine read_reference_model


!----------------------------------------------------------------------------------------!
! lib function linear_intp: lagrange interpolation
!----------------------------------------------------------------------------------------!
      subroutine linear_intp_dep(x1,y1,x2,y2,x,y)
      implicit none
      real(kind=4) ,intent(in):: x1,y1,x2,y2,x
      real (kind=4) ,intent(out) ::y
        y=y1+(y2-y1)/(x2-x1)*(x-x1)
       end subroutine linear_intp_dep
!----------------------------------------------------------------------------------------!
!  lib function locate.
!----------------------------------------------------------------------------------------!
     FUNCTION locate(xx,n,x)
!----------------------------------------------------------------------
! Given an array xx(1:N), and given a value x, returns a value j such
! that x is between xx(j) and xx(j+1). xx must be monotonic, either
! increasing or decreasing. j = 0 or j = N is returned to indicate
! that x is out of range.
!  NOTE:  the output j>=1, but j< n, so if  x=xx(n), the output j is still n-1, its a problem
!---------------------------------------------------------------------
       IMPLICIT NONE
       REAL(KIND=4), DIMENSION(n), INTENT(IN) :: xx
       REAL(KIND=4), INTENT(IN) :: x
       INTEGER(KIND=4) :: locate
       INTEGER(KIND=4) :: n,jl,jm,ju
       LOGICAL :: ascnd
       ascnd = (xx(n) >= xx(1))
        jl=0
        ju=n+1
       do
           if (ju-jl <= 1) exit
           jm=(ju+jl)/2
           if (ascnd .eqv. (x >= xx(jm))) then
              jl=jm
           else
              ju=jm
           end if
        end do
        if (x == xx(1)) then
          locate=1
        else if (x == xx(n)) then
           locate=n-1
        else
           locate=jl
        end if
       END FUNCTION locate

!----------------------------------------------------------------------------------------!
! function calculate_tref:  tref=ref period(secs) of model for dispersion correction.
!----------------------------------------------------------------------------------------!
      subroutine calculate_tref(kmaxRc,tRc,tref)
      IMPLICIT NONE
      integer kmaxRc
      real*8 tRc(kmaxRc)
      REAL*8 tref,tlog
       tlog=(log(tRc(1))+log(tRc(kmaxRc)) )/2
       tref=exp(tlog)

   END SUBROUTINE

