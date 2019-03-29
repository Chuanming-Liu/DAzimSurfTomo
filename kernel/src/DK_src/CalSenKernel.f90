! CODE FOR SURFACE WAVE  ANISOTROPY TOMOGRAPHY USING DISPERSION MEASUREMENT
! MAKE ANISOTROPIC KERNEL dc/dL  dc/dA
! VERSION: 1.0
! AUTHOR:
! CHUANMING LIU. chuanmingliu@foxmail.com
! HISTORY:
!        2016/10/22
       SUBROUTINE  CalSenKernel(nx,ny,nz,rmax, vel, kmaxRc, tRc,depz,&
       minthk,refmodel_file,ny_clc)
        IMPLICIT NONE
        integer nx,ny,nz,ny_clc
        real minthk
        real vel(nx,ny,nz)
        real vpz(nz),vsz(nz),rhoz(nz)
        integer ii,jj,k,i,nn,kk,j,z
        real depz(nz)
        integer kmaxRc     ! number of
        real*8 tRc(kmaxRc) ! periods array
        CHARACTER(*)refmodel_file
        integer mmax,rmax,rmax1
    	integer,parameter::NL=200
    	real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
        CHARACTER(LEN=50) key
        CHARACTER(LEN=6)::stri,strj,lonsrj
        integer nsublay(NL)

        REAL*8::LonSen_dcdL(nx,kmaxRc,rmax),LonSen_dcdA(nx,kmaxRc,rmax)
        REAL*8::LonSen_dcdC(nx,kmaxRc,rmax),LonSen_dcdF(nx,kmaxRc,rmax)
        REAL*8 :: dcR_dL(kmaxRc,rmax),dcR_dA(kmaxRc,rmax)
         REAL*8 :: dcR_dC(kmaxRc,rmax),dcR_dF(kmaxRc,rmax)
        INTEGER:: lay,senfid
        CHARACTER(LEN=100)senfile

! VARIABLE LIST
! kmaxRc: number of periods of specifciated types
! tRc: periods array
! depz: depth array
! minthk: num of refined layer
! ny_clc: cycle longtitude 1-ny
! nz: initial model grid num
! lay:inversion layer num
! jj=ny_clc: jj=1,ny
!---------------------------
        mmax=nz
        lay=nz-1
         jj=ny_clc

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
                call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax1,rdep,&
                rvp,rvs,rrho,rthk,nsublay)
                IF (rmax.NE.rmax1) STOP'Different rmax: Subroutine CalSenKernel.'

                 write(stri,'(I3.3)')ii
                 write(strj,'(I3.3)')jj
                  key='K'//TRIM(strj)//'_'//TRIM(stri)
                  write(*,'(a10)')key
                 call sensitivityMineos(NL,key,rmax,rdep,rvp,rvs,rrho,rthk,kmaxRc, tRc,&
                 refmodel_file,nsublay,nz, dcR_dL,dcR_dA,dcR_dC,dcR_dF)

                  LonSen_dcdL(ii,1:kmaxRc,1:rmax)=dcR_dL(1:kmaxRc,1:rmax)
                  LonSen_dcdA(ii,1:kmaxRc,1:rmax)=dcR_dA(1:kmaxRc,1:rmax)
                  LonSen_dcdC(ii,1:kmaxRc,1:rmax)=dcR_dC(1:kmaxRc,1:rmax)
                  LonSen_dcdF(ii,1:kmaxRc,1:rmax)=dcR_dF(1:kmaxRc,1:rmax)

            enddo
! OUTPUT
            write(lonsrj,'(I2.2)')ny_clc
            senfile='Lon_'//TRIM(lonsrj)//'_dcdL_dcdA.dat'
            senfid=98
            open(unit=senfid,file=senfile,action='write')
            write(senfid,*)ny_clc
            write(senfid,*)rmax
            write(senfid,*)(((LonSen_dcdL(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
            write(senfid,*)(((LonSen_dcdA(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
            write(senfid,*)(((LonSen_dcdC(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)
            write(senfid,*)(((LonSen_dcdF(j,i,z),z=1,rmax),i=1,kmaxRc),j=1,nx)

            close(98)
       END SUBROUTINE CalSenKernel

        subroutine sensitivityMineos(NL,key,rmax,rdep,rvp,rvs,rrho,rthk,&
                   kmaxRc, tRc,refmodel_file,nsublay,nz,dcR_dL,dcR_dA,dcR_dC,dcR_dF)
 !---------------------------------------------------------------------------------------!
!   The  code calls Mineos to calcualte the surface wave sensitivity kernel
!   by CH.M. Liu, 2016
!   NOTE:
!   1. use dcdA,dcdL in kernl_module to pass parameter, which is assigned in sub merge_kern of  buildG_dcdm1D_module_TI.f90
!   2. use buildG_dcdm1D_module to calculate sensitivity kernel with mineos
!----------------------------------------------------------------------------------------!
       implicit None
       CHARACTER(LEN=400)out_depth,outmdl,disp_file,infile,fout
       CHARACTER(*) key
       CHARACTER(LEN=400)fmerge
       INTEGER:: nz
       integer mmax,rmax,Ndep
       integer,INTENT(IN)::NL
       real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
        integer nsublay(NL)
       REAL*8 dcR_dL(kmaxRc,rmax),dcR_dA(kmaxRc,rmax)
       REAL*8 dcR_dC(kmaxRc,rmax),dcR_dF(kmaxRc,rmax)
       REAL:: Depth_array(NL)
       integer kmaxRc     ! number of
       real*8 tRc(kmaxRc) ! periods array

       REAL*8 dcdL(NL,NL),dcdA(NL,NL),dcdC(NL,NL),dcdF(NL,NL)
       CHARACTER(*)refmodel_file
       INTEGER:: i,j,k,z,jj
       real*4 :: mks=1000
       INTEGER:: ifid0=44,ifid_merge
       CHARACTER(LEN=12) :: fmtout='(1000E25.12)'
       CHARACTER(LEN=12):: fdepout='(1000F15.6)'
       real :: coe_vp(rmax),coe_rho(rmax),coe_A(rmax),coe_L(rmax),coe_C(rmax),coe_F(rmax)
       real*8 ::sen_Vs(kmaxRc,nz-1),sen_Gsc(kmaxRc,nz-1)
       real*8:: A_TIM(rmax),L_TIM(rmax)
!---------------------------------------------------------
! Main Variable List
! dz:  !  interpolation interval dz (km) in input model (interval of reference model is 2km in shallow ),  which depends on the period of data
!! INPUT:
!!         rmax: number of layers in the fined layered model, count number
!!         rdep:1. doesn't contain 0km; 2. stands for the lower interface depth;3. doesn't contain half space.
!!         rvp, rvs, rrho, rthk: the refined layered velocity model
!!         rmax, rvp, rvs, rrho, rthk : do not contain half space which is deleted in sub refinedGrid
!           kmaxRc: period num
!!         nsublay:[1:mmax-1] the refined layer num between every grid,mmax=nz
!!         nz: input initial model grid number
!!----------------------------
! Depth_array(1:1+rmax): mean the integral segment used by Mineos
! dcR_dL(kmaxRc,rmax)
!--------------------------------------------------------------------

! Step 1. Setting
! signal

!    output the model file of mineos format
       outmdl=TRIM(key)//'_model.card'
       disp_file=TRIM(key)//'_disp.dat'

! Step 2. Set mineos format model
       CALL mineos_model_output(refmodel_file,outmdl,&
       rmax,rdep,rvp,rvs,rrho,rthk,kmaxRc,tRc)

! output Integral depth segment
!  Note: dep_out contains 0km, but rdep does not contain 0Km

!       Open(Unit=42,File=out_depth,action='write')
!       write(Unit=42,FMT='(I3)')rmax !\ num of input depth segment.
!       write(Unit=42,FMT='(f9.1)')(rdep(i)*mks,i=1,rmax)
!       close(Unit=42)
       Depth_array(1)=0
       DO i=2,rmax+1
           Depth_array(i)=rdep(i-1)*mks
        END DO
! Step 3. the disp data of mineos format
        CALL mineos_disp_output(kmaxRc,tRc,disp_file)

! Step 4.   use Mineos to calculate the kernel
        WRITE (*,*)'USE buildG_dcdm1D with Mineos kernel to cal sensitivity kernel for A,L'
        CALL buildG_dcdm1D(key,disp_file,outmdl,Depth_array,rmax,dcdL,dcdA,dcdC,dcdF,NL)
!    load Mineos sensitivity kernel
!    For dcdL(i,j) i+,freq+, period- (sort in buildG_dcdm1D_module), depth ascending
!    dcR_dL: no sort,  we set the input data in period ascending order, depth ascending
        DO i=kmaxRc,1,-1
           dcR_dL(kmaxRc+1-i,1:rmax)=dcdL(i,1:rmax)
           dcR_dA(kmaxRc+1-i,1:rmax)=dcdA(i,1:rmax)
           dcR_dC(kmaxRc+1-i,1:rmax)=dcdC(i,1:rmax)
           dcR_dF(kmaxRc+1-i,1:rmax)=dcdF(i,1:rmax)
        ENDDO
! delete relevant files.

        OPEN(UNIT=ifid0,FILE=disp_file)
        CLOSE(UNIT=ifid0,STATUS='DELETE')

        OPEN(UNIT=ifid0,FILE=outmdl)
        CLOSE(UNIT=ifid0,STATUS='DELETE')

!  Step 5. output kernel for 1D model ii,jj
        fmerge=trim(key)//'_dL_dA_dC_dF.dat'
        ifid_merge=77
        WRITE(fmtout(2:5),'(I4.4)') rmax+1
        WRITE(fdepout(2:5),'(I4.4)')rmax+1

        OPEN(UNIT=ifid_merge,FILE=fmerge,ACTION='WRITE')
!  dc/dL absolute kernel
        WRITE(UNIT=ifid_merge,FMT=fdepout)0.0,(rdep(k),k=1,rmax)
        DO i=1,kmaxRc
            WRITE(ifid_merge,*) tRc(i),(dcR_dL(i,j),j=1,rmax)
        ENDDO
! dc/dA: absolute kernel
        WRITE(UNIT=ifid_merge,FMT=fdepout)0.0,(rdep(k),k=1,rmax)
        DO i=1,kmaxRc
            WRITE(ifid_merge,*)tRc(i),(dcR_dA(i,j),j=1,rmax)
        ENDDO
! dc/dC
        WRITE(UNIT=ifid_merge,FMT=fdepout)0.0,(rdep(k),k=1,rmax)
        DO i=1,kmaxRc
            WRITE(ifid_merge,*)tRc(i),(dcR_dC(i,j),j=1,rmax)
        ENDDO
! dc/dF
        WRITE(UNIT=ifid_merge,FMT=fdepout)0.0,(rdep(k),k=1,rmax)
        DO i=1,kmaxRc
            WRITE(ifid_merge,*)tRc(i),(dcR_dF(i,j),j=1,rmax)
        ENDDO

        CLOSE(UNIT=ifid_merge)


!-------------------------------------------------------------------------!
! Step 6. calculate relative kernel.

        END SUBROUTINE sensitivityMineos

         subroutine buildG_dcdm1D(key,fdata,fmodel,Depth_array,rmax,dcdL,dcdA,dcdC,dcdF,NL)
!-------------------------------------------------
! build sensitivity kernel at a given frequency
!
! H.Y. Yang, 2011
! Modified by CH.M. Liu, 2016
!------------------------------------------------
       USE data_type
       USE buildG_dcdm1D_module

       IMPLICIT NONE
       integer,intent(in)::NL
       INTEGER, PARAMETER :: itermax=20, nmax=6
       LOGICAL :: isderr  !\ if data contains data error
       INTEGER :: nfiles  !\ # of input data files
      CHARACTER(*),INTENT(in) :: fmodel, fdata
      CHARACTER(*),INTENT(IN)::key
      INTEGER:: rmax
      real:: Depth_array(NL)
      REAL*8 dcdL(NL,NL),dcdA(NL,NL)
      REAL*8 dcdC(NL,NL),dcdF(NL,NL)

       INTEGER :: njcom, jcoma(2), nbramax(2), nbra(nmax,2)
       REAL(KIND=8) :: freqa(2,2)
       LOGICAL :: isrunmode, ismulti, isddmks, isrelperturb

       INTEGER :: i,j,iter,nmodel_unit
       INTEGER:: ndata
! nmodel_unit: number of layer, rmax
!  set:
!# 'number of grids is power of 2? (T=yes, F=no)...>'  ismulti
!#  'Data is in kg-km-s(F) or mks(T)?  >'  isddmks
!# 'Relative perturbation kernel? (T)...>'  isrelperturb
!# 'Enter # of input file for data >' number of the input dispersion file >  nfiles
!# 'Enter filename of data...>'                    fdata
!# 'Enter initial guess model (*.card) ...> fmodel
!#  ''Enter depth segment file (in m)' ...>'  fdepth
! # 'Run mineos_kern to build table...(T=yes, F=no)>'  M_000_Sph.fre>
!#  # 'Enter the xth iteration...>'----Flag 000 > iter
       isderr=.TRUE.
       ismulti=.FALSE.
       isddmks=.FALSE.
       isrelperturb=.FALSE.
       nfiles=1

      isrunmode=.TRUE.  !\ 'Run mineos_kern to build table...>'
       iter=0 !\xth iteration
!%-----------------------------------------------------------------------------------------%!
!  process
!%-----------------------------------------------------------------------------------------%!
!       ndata=0 !/ must set ndata=0 here to avoid compiler-dependent error
        ndata=0
       CALL read_data(fdata,isddmks,ndata)
       CALL sortdata(ndata)
       CALL datarange(njcom, jcoma, freqa, nbramax, nbra,ndata)
!  WRITE(*,*) 'nbra=', ((nbra(i,j),i=1,nbramax(j)+1),j=1,njcom)

!      Enter depth segment file (in m)' Depth_array
       write(*,*)'Please ignore: forrtl: warning (406)!'
       nmodel_unit=rmax
       CALL get_dhat_kern(key,fmodel,iter,njcom,jcoma,freqa,nbramax,nbra,&
                    isrunmode,ismulti,NL,dcdL,dcdA,dcdC,dcdF,Depth_array,nmodel_unit)

END SUBROUTINE buildG_dcdm1D



        subroutine mineos_disp_output(kmaxRc,tRc,disp_file)
!-------------------------------------------------
! Convert the  input disp data to Mineos input disp file
!
! CH.M. Liu, 2016
!-----------------------------------------------------------------------------!
!                         mineos  data format                                               !
! format: /mode_type /data_type/  branch(n)  / freq(Hz) /  phase v/ error/
!   CPS:    /indexRL/ indexCG /  indexMode- 1/period(i)/0/0/
!                   CPS-----> Mineos
!   indexRL=1(Love), 2(Rayleigh)-> mode_type= 3(Rayleigh),2(Love)
!   indexCG=0(phase), 1(group)->data_type= 1(phase )
!  indexMode=1, fundamental mode-> branch(n) =0 (fundamental mode); 1(1st overtone)  etc.
!-----------------------------------------------------------------------------!
        implicit none
        character(Len=*),intent(in)::  disp_file
        integer,intent(in):: kmaxRc     ! number of period
        real*8,intent(in):: tRc(kmaxRc) ! periods array
        integer:: i,j
        real:: phaseV,Verror
        INTEGER:: mode_type,data_type,branch
        mode_type=3
        data_type=1
        branch=0
        phaseV=0.0
        Verror=0.0
          open(unit=23,file=disp_file,action='write')
          write(unit=23,FMT='(I4)')kmaxRc
          write(unit=23,FMT='(I3,I3,I3,F9.5 ,F7.3,F7.3)')(mode_type,data_type,branch,1/tRc(i),phaseV,Verror,i=1,kmaxRc)
          close(unit=23)
         end subroutine mineos_disp_output


       subroutine refineGrid2LayerMdl(minthk0,mmax,dep,vp,vs,rho,&
                  rmax,rdep,rvp,rvs,rrho,rthk,nsublay)
!!--------------------------------------------------------------------c
!!refine grid based model to layerd based model
!!INPUT:   minthk: is the minimum thickness of the refined layered model
!                   minthk0: is the number of cut time, minthk0=1, 2 layers
!!         mmax: number of depth grid points in the model
!!         dep, vp, vs, rho: the depth-grid model parameters
!!OUTPUT:
!!         rmax: number of layers in the fined layered model, count number
!!         rdep:1. doesn't contain 0km; 2. stands for the lower interface depth
!!         rvp, rvs, rrho, rthk: the refined layered velocity model
!!         nsublay:[1:mmax-1] the refined layer num between every grid
!! NOTE: need to consider whether need to keep the effect of halp space
!! Note: dep contains 0km, but rdep does not contain 0Km
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
! commented by CH.M.L, Oct,2016
!  Considering the  refined model is  midpoint of the layer point, it don't need half space.
!        k = k + 1
!        rthk(k) = 0.0
!        rvp(k) = vp(mmax)
!        rvs(k) = vs(mmax)
!        rrho(k) = rho(mmax)
!	rdep(k) = dep(mmax)

        rmax = k

!!       do i = 1, mmax
!!          write(*,*) dep(i),vp(i),vs(i),rho(i)
!!       enddo
!!       print *, '---------------------------------'
!!       do i = 1, rmax
!!          write(*,*) rdep(i),rthk(i),rvp(i),rvs(i),rrho(i)
!!       enddo

        return
        end
