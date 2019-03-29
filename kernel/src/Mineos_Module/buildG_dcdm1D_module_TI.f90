MODULE buildG_dcdm1D_module
!------------------------------------
! H.Y. Yang, Mar 2015
!  debug opening file of fmerge
!
! Developed by H.Y. Yang, 2010
!  -calculate sensitivity kernel of phase velocity
!
! execute external program mineos_kern (phase velocity kernel)
!------------------------------------
  USE data_type
  USE math_par_8b, ONLY : pi_r8
  USE sort_module
  USE interpol_module, ONLY : polint, locate
  USE modenorm_module  !/ rhon, rn, gn, bigg, wn2

  IMPLICIT NONE
  CHARACTER(LEN=500) :: mineos_kern_path ='/home/lance/lance/project/DirecInv_kernel/source_new/Mineos_kern/&
  mineos_kern_dcdL_TI'
  LOGICAL, PARAMETER :: isalpha=.TRUE., isbeta=.TRUE., isrho=.TRUE. !\ for output
  INTEGER, PARAMETER :: iflag_phs=1, iflag_grp=2, iflag_sph=3, iflag_tor=2
  INTEGER, PARAMETER :: nrmax=50000, nhead=6,nkern=7
  REAL(KIND=r8b), PARAMETER :: km2m=1000.0_r8b, twopi=2.0_r8b*pi_r8
  REAL(KIND=8),PARAMETER :: m2km=1.d-3
!  kernel type: 1:dw/dv, 2: dc/dv 3. dc/dL
  INTEGER, PARAMETER :: iflag_kernel=3

  !%=============================%
  !|  Subroutine read_data       |
  !%=============================%
  INTEGER, PARAMETER :: ndatamax=500
  TYPE datafmt
    INTEGER ::        jcom  !/mode type =2:tor, =3:sph
    INTEGER ::        icu   !/phase or group velocity =1:phs, =2, grp
    INTEGER ::        nbran !/n branch
    REAL(KIND=r8b) :: freq  !/frequency in Hz
    REAL(KIND=r8b) :: dd    !/data in m/s
    REAL(KIND=r8b) :: derr  !/data error in m/s
  END TYPE datafmt
  TYPE(datafmt) ,ALLOCATABLE:: data0(:)
! COMMENTED BY C.M. L
!  INTEGER :: ndata=0
  !%=============================%
  !| Subroutine open_fmerge_gram |
  !%=============================%
  TYPE minfotype
  !/ Note that xb to xe is in increasing depth
    REAL(KIND=r8b) :: xb
    REAL(KIND=r8b) :: xe
    REAL(KIND=r8b) :: dx !\grid spacing for model to be inverted in m
  END TYPE minfotype
  TYPE modeltype
    REAL(KIND=r8b) :: rr
    REAL(KIND=r8b) :: rho
    REAL(KIND=r8b) :: vp
    REAL(KIND=r8b) :: vs
    REAL(KIND=r8b) :: Qmu
    REAL(KIND=r8b) :: Qa
  END TYPE modeltype


  TYPE(minfotype) :: minfo
  !/ nmodel_unit=NINT((xe-xb)/dx)+1
  TYPE(modeltype), ALLOCATABLE :: model_1d( : ) !/in mks
  INTEGER :: nmodel
  REAL(KIND=r8b), ALLOCATABLE :: gram_row( : ) !\gram_row(nmodel)
  INTEGER, ALLOCATABLE ::  irbds( : , : )      !\irbds(2,nmodel_unit)

! Added by C.M L  default rayleigh wave,fundamental mode
    !%=============================%
  !|  Subroutine load_integral_depth       |
  !%=============================%
  TYPE  depthtype
   REAL(KIND=r8b)::  up !/ depth segment of upper boundary
   REAL(KIND=r8b):: dw !/ depth segment of downer boundary
  ENDTYPE depthtype

CONTAINS
  SUBROUTINE read_data(fin,isunit_mks,ndata)
  !--> data0, ndata
  ! use km2m
  !-----------------------------------------------
  ! read surface wave data and change unit to mks if necessary
  ! L1   # of measurement
  ! L?   jcom, icu, nbran, freq(Hz), dd(m/s), derr(m/s)
  !-----------------------------------------------
  CHARACTER(LEN=128), INTENT(IN) :: fin
  LOGICAL, INTENT(IN) :: isunit_mks
! Added by C.M. L
  INTEGER,INTENT(OUT) :: ndata
  INTEGER, PARAMETER :: fid=20
  INTEGER :: nd
  INTEGER :: Ier, i, ib, ie
  LOGICAL :: isofr
! Added by C.M. L may cause   repeat Data.
  ndata=0
  IF (ALLOCATED(data0)) DEALLOCATE(data0)


  OPEN(UNIT=fid,FILE=fin,ACTION='READ',STATUS='OLD',IOSTAT=Ier)
  IF (Ier>0) STOP 'Reading data error, Terminated!'
  READ(UNIT=fid,FMT=*) nd
  ib=ndata+1
  ndata=ndata+nd

  ALLOCATE(data0(nd),STAT=Ier)
  IF (Ier/=0) STOP 'ALLOCATE error of model_1d in sub read_data'

  IF (ndata> ndatamax) STOP 'ndatamax must be larger than ndata, Terminated!'
  ie=ndata
  READ(UNIT=fid,FMT=*) (data0(i),i=ib,ie)
  CLOSE(UNIT=fid)
  isofr = &
  ANY((data0(ib:ie)%jcom/=iflag_tor .AND. data0(ib:ie)%jcom/=iflag_sph) .OR. &
      (data0(ib:ie)%icu/=iflag_phs  .AND. data0(ib:ie)%icu/=iflag_grp))
  IF (isofr) STOP 'jcom and icu out of range, Terminated!'
  IF (.NOT.isunit_mks) THEN
    data0(ib:ie)%dd=data0(ib:ie)%dd*km2m
    data0(ib:ie)%derr=data0(ib:ie)%derr*km2m
  END IF
!   OPEN(UNIT=fid,FILE=fin)
!   CLOSE(UNIT=fid,STATUS='DELETE')
  END SUBROUTINE

  SUBROUTINE sortdata(ndata)
  ! -> data0
  ! use global variables: data0 and ndata from sub read_data
  ! USE sort_module
  !---------------------------------------------------------
  ! sort data0 according to jcom, icu, nbran, freq in turn
  !---------------------------------------------------------
  ! Added by C.M. L
  INTEGER,INTENT(INOUT) :: ndata

  INTEGER, ALLOCATABLE :: index_array( : ), jtype( : )
  REAL, ALLOCATABLE :: rtype( : )
  REAL :: scalen
  INTEGER :: jseg(4), nseg
  INTEGER :: Ier, i, i1, j, ib, ie

  !%--------------------------------------------------------%
  !| sort according to jtype which is jcom multiplying icu  |
  !%--------------------------------------------------------%
  ALLOCATE(jtype(ndata),index_array(ndata),STAT=Ier)
  IF (Ier>0) STOP 'Allocate error of index_array in sub sortdata, STOP!'
  jtype = data0(1:ndata)%jcom*data0(1:ndata)%icu
  CALL indexx(jtype,index_array)
  jtype=jtype(index_array)

  !%-------------------------------------------------------%
  !| find segment of jtype -> nseg, jseg                   |
  !%-------------------------------------------------------%
  nseg=1
  jseg=0
  j=1
  jseg(j)=1
  i=1
  DO
    i1=i+1
    IF (i1>ndata) EXIT
    IF (jtype(i1)-jtype(i) /= 0) THEN
      j=j+1
      IF (j > 4) STOP &
       'rank of jseg in sub sortdata should be greater &
        than current value, Terminated!'
       jseg(j)=i1
       nseg=nseg+1
      END IF
    i=i1
  END DO
  DEALLOCATE(jtype)

  !%-----------------------------------------------------------------%
  !| Within each segment, arrange data0 based on frequency+          |
  !|   upperbound(freq)*nbran                                        |
  !%-----------------------------------------------------------------%
  data0(1:ndata)=data0(index_array)
  ALLOCATE(rtype(ndata),STAT=Ier)
  IF (Ier > 0) STOP 'allocate error of rtype in sub sortdata'
  scalen=CEILING(MAXVAL(data0(1:ndata)%freq))*1.0
  rtype=scalen*data0(1:ndata)%nbran+REAL(data0(1:ndata)%freq)

  DO i = 1, nseg
    ib=jseg(i)
    IF ( i == nseg) THEN
      ie=ndata
    ELSE
      ie=jseg(i+1)-1
    END IF
    CALL indexx(rtype(ib:ie), index_array(ib:ie))
    index_array(ib:ie)=index_array(ib:ie)+ib-1
  END DO
  data0(1:ndata)=data0(index_array)
  IF (ANY(data0(1:ndata)%jcom-data0(1)%jcom /= 0) ) THEN
    CALL indexx(data0(1:ndata)%jcom, index_array)
    data0(1:ndata)=data0(index_array)
  END IF

!  WRITE(*,FMT='(I2,I2,I2,3F9.3)') (data0(i),i=1,ndata)

  DEALLOCATE(index_array,rtype)
  END SUBROUTINE sortdata

  SUBROUTINE datarange(nseg,jcoma,freqa,nbramax,nbra,ndata)
  ! use data0, ndata
  !---------------------------------------------------------------
  ! find frequency range (freqa) and maximal branch (nbramax), and
  !  # of modes for a given specific branch
  ! for corresponding mode type (jcoma)
  !---------------------------------------------------------------
  INTEGER, INTENT(OUT) :: nseg, jcoma(2), nbramax(2)
  REAL(KIND=r8b), INTENT(OUT) :: freqa(2,2)
  INTEGER, INTENT(OUT) :: nbra(:,:)
  ! Added by C.M. L
  INTEGER,INTENT(IN) :: ndata

  INTEGER :: nhit, ndata1, idx, jsegb(2), jsege(2), i, k, nn

  IF (SIZE(nbra,DIM=1)<MAXVAL(data0(1:ndata)%nbran)+1) &
    STOP 'Dimension of nbra is too small, Terminated!'
  nbra=0
  jcoma=0
  freqa=0.0_r8b
  ndata1=ndata-1
  nhit=COUNT(data0(2:ndata)%jcom-data0(1:ndata1)%jcom==0)
  IF (nhit == ndata1) THEN
    nseg = 1
    jcoma(1)= data0(1)%jcom
    freqa(1,1)=MINVAL(data0(1:ndata)%freq)
    freqa(2,1)=MAXVAL(data0(1:ndata)%freq)
    nbramax(1)=MAXVAL(data0(1:ndata)%nbran)
    DO k = 1, nbramax(1)+1
      nn=k-1
      nbra(k,1)=COUNT(data0(1:ndata)%nbran==nn)
    END DO
  ELSE IF (nhit == ndata1-1) THEN
    idx=MINLOC(data0(1:ndata1)%jcom, &
               MASK=data0(2:ndata)%jcom-data0(1:ndata1)%jcom==1,&
               DIM=1)+1
    jsegb=(/ 1, idx /)
    jsege=(/ idx-1, ndata /)
    nseg = 2
    DO i = 1, nseg
      jcoma(i)=data0(jsegb(i))%jcom
      freqa(1,i)=MINVAL(data0(jsegb(i):jsege(i))%freq)
      freqa(2,i)=MAXVAL(data0(jsegb(i):jsege(i))%freq)
      nbramax(i)=MAXVAL(data0(jsegb(i):jsege(i))%nbran)
      DO k = 1, nbramax(i)+1
        nn=k-1
        nbra(k,i)=COUNT(data0(jsegb(i):jsege(i))%nbran==nn)
      END DO
    END DO
  ELSE
    STOP 'data are not sorted as we wish, Terminate!'
  END IF
  IF (SUM(nbra)/=ndata) STOP 'error in calculating nbra in sub datarange!'
  END SUBROUTINE datarange

  SUBROUTINE get_dhat_kern(key,fmodel, iter, njcom, jcoma, freqa, nbramax, nbra,&
                           isxmineos, ismulti,NNL,dcdL,dcdA,dcdC,dcdF,Depth_array,nmodel_unit)
  !  delete is_relperturb in input parameter list
  ! call intrinsic function "system" to run mineos_kern
  ! call split_branch
  ! use km2m
  ! USE interpol_module, ONLY : polint, locate
  CHARACTER(LEN=300), INTENT(IN) :: fmodel
  INTEGER, INTENT(IN) :: iter
  INTEGER, INTENT(IN) :: njcom, jcoma( : ), nbramax( : ), nbra( : , : )
  LOGICAL, INTENT(IN) :: isxmineos, ismulti
  CHARACTER(*),INTENT(IN)::key
  LOGICAL:: is_relperturb !\ deleted
  REAL(KIND=r8b), INTENT(IN) :: freqa(:,:)  !\frequency range in Hz
! Added by L.C.M
  INTEGER,INTENT(IN):: NNL
  REAL(KIND=r8b) dcdA(NNL,NNL),dcdL(NNL,NNL)
  REAL(KIND=r8b) dcdC(NNL,NNL),dcdF(NNL,NNL)
  INTEGER,INTENT(IN):: nmodel_unit
  REAL(KIND=r4b),INTENT(IN)::Depth_array(NNL)

  REAL(KIND=r8b), PARAMETER :: eps=1E-15_r8b, grav=1E+03_r8b, &
      lmin=1.0_r8b, lmax=20000.0_r8b, dl=0.5_r8b, &
      dfmhz=1.0_r8b
  INTEGER, PARAMETER :: ifid0=39, ifid_merge=41, ifid_txt=42
  CHARACTER(LEN=4) :: id='_000', mtype
  CHARACTER(LEN=600) :: finp, frename, funname, fkern, ftmp, string, &
                        fout, fsplit, ftmp1, fmerge
  CHARACTER(LEN=4) :: cnb='_n00'
  REAL(KIND=r8b),DIMENSION(SIZE(freqa,DIM=1),SIZE(freqa,DIM=2)) :: &
       freqa_mhz
  LOGICAL :: isexist, isphs, isnewjcom, isfalse
  REAL(KIND=r8b), ALLOCATABLE :: ltbl( : ), ftbl( : )
  REAL(KIND=r8b) :: f_target,l_target,dl_target, r8btmp3(3), truc_depm
  INTEGER :: nl, nmodes(SIZE(nbra,DIM=1))
  INTEGER :: i,j,k, Ier, idata, idatab, idatae, nbra_sum, ibrk, ibrk1, ibrk2
  INTEGER:: nmin,nmax

   TYPE(depthtype):: depth_seg(NNL)
! Load depth array for integral
  depth_seg(1:nmodel_unit)%up=Depth_array(1:nmodel_unit)
  depth_seg(1:nmodel_unit)%dw=Depth_array(2:nmodel_unit+1)

! for we only  process fundamental modes, so set nbranch =0,nmin,nmax=0
  nmin=0
  nmax=0

!  truc_depm=truc_dep()*km2m
   truc_depm=depth_seg(nmodel_unit)%dw

  WRITE(*,*) 'truc_dep in truc_dep(km)=', truc_depm/km2m
  nl=CEILING((lmax-lmin)/dl)
  ALLOCATE(ltbl(nl),ftbl(nl),STAT=Ier)
  IF (Ier/=0) STOP 'Allocate error in ltbl in sub get_dhat_kern, Terminate!'

  freqa_mhz=freqa*1000.0_r8b !\ frequency range in mhz
  freqa_mhz(1,:)=freqa_mhz(1,:)-50.0_r8b
  freqa_mhz(2,:)=freqa_mhz(2,:)+5.0_r8b
  DO i = 1, njcom
    IF (freqa_mhz(1,i) <= 0.0_r8b) freqa_mhz(1,i)= 1.0_r8b
  END DO
  WRITE(id(2:4),'(I3.3)') iter

  fmerge=TRIM(key)//'Gramdd'//TRIM(id)//'.Gd'
  OPEN(UNIT=ifid_txt,FILE=TRIM(fmerge)//'used')

  nbra_sum=0
  DO i = 1, njcom
   isnewjcom = .TRUE.
   !%-----------------------------------%
   !| create input file for mineos_kern |
   !%-----------------------------------%
    IF (jcoma(i) == iflag_sph) THEN
      mtype='_Sph'
    ELSE IF (jcoma(i) == iflag_tor) THEN
      mtype='_Tor'
      STOP 'Love wave-Tor type input.'
    END IF
    ftmp=TRIM(key)//'M'//TRIM(id)//TRIM(mtype)
    finp=TRIM(ftmp)//'.inp'
    frename=TRIM(ftmp)//'.fre' ! modal table
    funname=TRIM(ftmp)//'.fun'
    fkern=TRIM(ftmp)//'.kern'
    fout=TRIM(ftmp)//'.out'

    OPEN(UNIT=ifid0,FILE=finp,ACTION='WRITE')
      WRITE(UNIT=ifid0,FMT='((A)/(A)/(A)/2(1X,E22.14)/I3/,&
         E13.6,1X,E13.6,2(1X,F8.3),I5,I5, 1X, E13.6/,&
         (A)/,I4)') &
      TRIM(fmodel), TRIM(frename), TRIM(funname),&
      eps, grav, &
      jcoma(i),&
      lmin, lmax, freqa_mhz(1,i), freqa_mhz(2,i),nmin,nmax, dl, &
      TRIM(fkern), iflag_kernel    !\iflag_kernel=3,dc/dL
    CLOSE(UNIT=ifid0)

    string=TRIM(mineos_kern_path)//' < '//TRIM(finp)//' >'//TRIM(fout)

    !%-------------------------------------%
    !| excute external program mineos_kern |
    !%-------------------------------------%
    IF ( isxmineos ) THEN
       Write(*,*)'Call mineos_kernel to make modal table.'
      CALL system(string)
      CALL system('sleep 1')
      WRITE(*,'(A)') 'Finishing mineos_kern for '//TRIM(ftmp)
    END IF

    !%------------------------------------------------------%
    !| split mineos output file according to their branches | make M_000_Sphn000.split
    !%------------------------------------------------------%
    INQUIRE(FILE=frename,EXIST=isexist)
    IF ( isexist ) CALL split_branch(nbramax(i),ftmp,nmodes)
    WRITE(*,'(A)') 'Finishing split_branch for '//TRIM(ftmp)
    !%-------------------------------------------------------%!
    ! | delete the M_000_Sph.* file for too big, especially for M_000_Sph.kern file

!    OPEN(UNIT=ifid0,FILE=frename)
!    CLOSE(UNIT=ifid0,STATUS='DELETE')!\ M_000_Sph.fre

    OPEN(UNIT=ifid0,FILE=funname)!\ M_000_Sph.fun
    CLOSE(UNIT=ifid0,STATUS='DELETE')
    OPEN(UNIT=ifid0,FILE=fkern) !\ M_000_Sph.kern
    CLOSE(UNIT=ifid0,STATUS='DELETE')
    OPEN(UNIT=ifid0,FILE=fout)
    CLOSE(UNIT=ifid0,STATUS='DELETE') !\ M_000_Sph.out
    OPEN(UNIT=ifid0,FILE=finp)
    CLOSE(UNIT=ifid0,STATUS='DELETE') !\ M_000_Sph.in
    !%------------------------------------------------------%
    !| polyfit each branch, i.e. ltbl as a function of ftbl |
    !%------------------------------------------------------%
    DO k = 1, nbramax(i)+1
      IF (k.GT.1) STOP'Higher branch input.'
      IF (nbra(k,i) > 0 ) THEN
        idatab=nbra_sum+1
        nbra_sum=nbra_sum+nbra(k,i)
        idatae=nbra_sum
        WRITE(cnb(3:4),'(I2.2)') k-1
        fsplit=TRIM(ftmp)//cnb//'.split'
        write(*,*)'SPLIT:',fsplit
        OPEN(UNIT=ifid0,FILE=fsplit,ACTION='READ',STATUS='OLD')
          READ(UNIT=ifid0,FMT=*) &
           (ltbl(j), ftbl(j), r8btmp3, j=1,nmodes(k)) !/ya,xa (in mhz)

!        CLOSE(UNIT=ifid0,STATUS='DELETE') !\delete M_000_Sph_n00.split
        CLOSE(UNIT=ifid0) !\delete M_000_Sph_n00.split
        ftbl=ftbl/1000.0_r8b  !/ in Hz

        DO idata = idatab, idatae
          isphs = .FALSE.
          IF ( data0(idata)%icu == iflag_phs ) isphs = .TRUE.
          f_target=data0(idata)%freq
          ibrk=locate(ftbl(1:nmodes(k)),f_target)
          ibrk1=ibrk-1
          ibrk2=ibrk+2
          CALL polint(ftbl(ibrk1:ibrk2),ltbl(ibrk1:ibrk2), &
                      f_target,l_target,dl_target)
!          WRITE(*,'(3F14.6)') l_target, f_target*1000.0_r8b, dl_target

         !%-----------------------------------------------%
         !| create mineos_kern input file and then run    |
         !%-----------------------------------------------%
          f_target=f_target*1000.0_r8b !/ in mhz
          ftmp1=TRIM(key)//'kern_recal'
          finp=TRIM(ftmp1)//'.inp'
          frename=TRIM(ftmp1)//'.fre'
          funname=TRIM(ftmp1)//'.fun'
          fkern=TRIM(ftmp1)//'.kern'
          fout=TRIM(ftmp1)//'.out'
          OPEN(UNIT=ifid0,FILE=finp,ACTION='WRITE')
            WRITE(UNIT=ifid0,FMT='((A)/(A)/(A)/,2(1X,E22.14)/,I3/,&
           E13.6,1X,E13.6,2(1X,F8.3),I5,I5, 1X, E13.6/,&
           (A)/,I4)') &
             TRIM(fmodel), TRIM(frename), TRIM(funname),&
            eps, grav, &
             jcoma(i),&
             l_target, l_target, f_target-dfmhz, f_target+dfmhz,nmin,nmax,dl, &
             TRIM(fkern), iflag_kernel
          CLOSE(UNIT=ifid0)
          string=TRIM(mineos_kern_path)//' < '//TRIM(finp)//' >'//TRIM(fout)
          CALL system(string)
          IF ( isnewjcom ) THEN

            CALL open_fmerge_gram(ismulti,ifid_merge,fmerge,fkern,truc_depm,NNL,depth_seg,nmodel_unit)
            isnewjcom = .FALSE.
          END IF
          CALL merge_kern(frename,fkern,ifid_merge,ifid_txt,&
                          f_target,isphs,idata,isfalse,NNL,dcdL,dcdA,dcdC,dcdF,nmodel_unit)
          IF ( isfalse ) &
            WRITE(*,'(A,I5)') 'Fail in writing kernel for data: ', idata
          !%----------------------------%
          !| delete output from mineos  |
          !%----------------------------%
          OPEN(UNIT=ifid0,FILE=fkern) ! kern_recal.kern
          CLOSE(UNIT=ifid0,STATUS='DELETE')
          OPEN(UNIT=ifid0,FILE=frename) ! kern_recal.fre
          CLOSE(UNIT=ifid0,STATUS='DELETE')
          OPEN(UNIT=ifid0,FILE=funname) ! kern_recal.fun
          CLOSE(UNIT=ifid0,STATUS='DELETE')
          OPEN(UNIT=ifid0,FILE=fout) ! kern_recal.out
          CLOSE(UNIT=ifid0,STATUS='DELETE')
          OPEN(UNIT=ifid0,FILE=finp) ! kern_recal.inp
          CLOSE(UNIT=ifid0,STATUS='DELETE')



        END DO !/ end of idata
!      CALL kernel_output(ifid_merge,fmerge)

      END IF!/ if (nbra(k,i) > 0 )

    END DO  !/ end of branch , cycle on k
    ! delete Modal table
    OPEN(UNIT=ifid0,FILE=frename)!\ M_000_Sph.fre
    CLOSE(UNIT=ifid0,STATUS='DELETE')

  END DO !/ end of jcom, cycle on i






  CLOSE(UNIT=ifid_txt) ! Gramdd_000.Gdused

  IF (ALLOCATED(data0)) DEALLOCATE(data0)
  DEALLOCATE(model_1d)
  DEALLOCATE(ltbl,ftbl)
  END SUBROUTINE get_dhat_kern

  SUBROUTINE split_branch(nbran,fname,nmodes)
  ! split modes according to their branches
  ! only read ASCII file (*.fre) from mineos
  INTEGER, INTENT(IN) :: nbran          !\ maximal radial order
  CHARACTER(LEN=*), INTENT(IN) :: fname !\ output filename from mineos
  INTEGER, INTENT(OUT) :: nmodes( : )   !\ modes number for every radial order

!  REAL(KIND=r8b), PARAMETER :: tol=1E-5_r8b
  REAL(KIND=r8b), PARAMETER :: tol=1E-4_r8b
  INTEGER :: nn
  REAL(KIND=r8b) :: ll, wrad, fmhz, tt, gcom, qmod, wdiff
  INTEGER :: nsplit
  INTEGER, ALLOCATABLE :: ncount(:), njump( : )
  REAL(KIND=r8b), ALLOCATABLE :: ll0( : )
  CHARACTER(LEN=128) :: frename, cjcom*2
  CHARACTER(LEN=128), ALLOCATABLE :: fresplit( : )
  CHARACTER(LEN=4) :: cnb='_n00'
  INTEGER :: Ier, fid0, fid1, i
  fid0=21
  frename=TRIM(fname)//'.fre'
  nsplit=nbran+1

  !%-----------------------------%
  !| open split files            |
  !%-----------------------------%
  ALLOCATE(fresplit(nsplit),ncount(nsplit),ll0(nsplit),njump(nsplit),STAT=Ier)
  IF (Ier/=0) STOP 'allocate error of fresplit in sub split_branch'
  DO i = 1, nsplit
    fid1=fid0+i
    WRITE(cnb(3:4),'(I2.2)') i-1
    fresplit(i)=TRIM(fname)//cnb//'.split'
    OPEN(UNIT=fid1,FILE=fresplit(i),ACTION='WRITE')
  END DO

  !%----------------------------------------------%
  !| split modes according to their radial number |
  !%----------------------------------------------%
  ll0=0
  ncount=0
  njump=0
  OPEN(UNIT=fid0,FILE=frename,ACTION='READ')
    DO
      READ(UNIT=fid0,FMT=*,IOSTAT=Ier) &
         nn, cjcom, ll, wrad, fmhz, tt, gcom, qmod, wdiff
      IF (Ier>0) THEN
        STOP 'Reading mineos output error in sub split_branch!'
      ELSE IF (Ier < 0 ) THEN
        EXIT
      END IF
      IF ((nn<=nbran) .AND. (ABS(wdiff)<tol)) THEN
        i=1+nn
        fid1=fid0+i
        ncount(i)=ncount(i)+1
        WRITE(UNIT=fid1,FMT='(E14.6,4(1X,E16.7))') ll, fmhz, gcom, qmod, wdiff
        IF ((ll-ll0(i) > 1.0_r8b) .AND. (ncount(i) > 1)) njump(i)=njump(i)+1
        ll0(i)=ll
      END IF
    END DO
  CLOSE(UNIT=fid0)

  !%-------------------------------------------%
  !| close split files and delete empty files  |
  !%-------------------------------------------%
  DO i = 1, nsplit
    fid1=fid0+i
    IF (ncount(i) == 0) THEN
      CLOSE(UNIT=fid1,STATUS='DELETE')
    ELSE
      IF (njump(i)/=0) WRITE(*,'(A,A,A,I3,A)') &
        '   Incomplete mode table for ', TRIM(fresplit(i)),&
        ' ...# of jumping= ', njump(i)
      CLOSE(UNIT=fid1)
    END IF
  END DO
  nmodes=0
  nmodes(1:nsplit)=ncount

  DEALLOCATE(fresplit,ncount)
  END SUBROUTINE split_branch

  SUBROUTINE open_fmerge_gram(ismulti,ifid_merge,fmerge,fkern,truc_depm,NNL,depth_seg,nmodel_unit)
  ! NOTE: deleted is_relative_perturb in parameter list.
  !-> nmodel, nmodel_unit, gram_row(blank), minfo, irbds, model_1d
  ! USE interpol_module, ONLY: locate
  ! USE USE modenorm_module, ONLY : rn
  ! use global variables: isbeta, isalpha, isrho, nrmax,
  !-----------------------------------------------------------
  ! open fmerge at the first time.
  ! calculate # of inversed model grids -> nmodel_unit, nmodel
  ! and then allocate gram_row for accessing directly thereafter.
  !-----------------------------------------------------------
  LOGICAL, INTENT(IN) :: ismulti
  INTEGER, INTENT(IN) :: ifid_merge
  CHARACTER(LEN=*), INTENT(IN) :: fkern, fmerge
  REAL(KIND=r8b), INTENT(INOUT) :: truc_depm
!  LOGICAL, INTENT(IN) :: is_relative_perturb
! Added by C.M. L
  INTEGER, INTENT(IN)  :: NNL
   TYPE(depthtype), INTENT(IN) :: depth_seg(NNL)
   INTEGER,INTENT(IN) :: nmodel_unit

  INTEGER, PARAMETER :: fid1=21, ndismax=20
  REAL(KIND=r8b) :: zero_app = EPSILON(1.0_r8b)*10.0_r8b
  CHARACTER(LEN=128) :: fmergepp
  INTEGER :: iolenGde, nlevel, npt
  REAL(KIND=r8b) :: xbound(2), rbound(2), cdep,rbound_up,rbound_dw
  REAL :: tdum(nrmax,9)
  INTEGER :: n, nic, noc
  INTEGER :: i, j, Ier
  INTEGER :: itmp
  REAL :: rtmp
  LOGICAL :: isexist
  !%------------------------------------------------------------%
  !| 1. Read the header of kernel file for info on reference model |
  !| (the following variables are not normalized)               |
  !%------------------------------------------------------------%
  OPEN(UNIT=fid1,FILE=fkern,STATUS='OLD',ACTION='READ',&
      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
!    READ(UNIT=fid1) jcom, wmin, wmax, flmin, flmax, wgrav
!    READ(UNIT=fid1) n,nic,noc,ifanis,trefl,((tdum(i,j),i=1,n),j=1,9)
    READ(UNIT=fid1) itmp, rtmp, rtmp, rtmp, rtmp, rtmp
    READ(UNIT=fid1) n, nic, noc, itmp, rtmp, ((tdum(i,j),i=1,n),j=1,9)
  CLOSE(UNIT=fid1)
  IF (ALLOCATED(model_1d)) DEALLOCATE(model_1d)
  ALLOCATE(model_1d(n),STAT=Ier)
  IF (Ier/=0) STOP 'ALLOCATE error of model_1d in sub open_fmerge_gram'
  DO i = 1, n
    model_1d(i)%rr =REAL(tdum(i,1),KIND=r8b)
    model_1d(i)%rho=REAL(tdum(i,2),KIND=r8b)
    model_1d(i)%vp =REAL(tdum(i,3),KIND=r8b)
    model_1d(i)%vs =REAL(tdum(i,4),KIND=r8b)
    model_1d(i)%Qa =REAL(tdum(i,7),KIND=r8b) ! here Qa is Qkappa
    model_1d(i)%Qmu=REAL(tdum(i,8),KIND=r8b)
    rtmp=4.0_r8b/3.0_r8b*(model_1d(i)%vs/model_1d(i)%vp)**2
    IF (model_1d(i)%Qmu > 0.0_r8b) THEN
      model_1d(i)%Qa =(1.0_r8b-rtmp)/model_1d(i)%Qa+rtmp/model_1d(i)%Qmu
      model_1d(i)%Qa = 1.0_r8b/model_1d(i)%Qa
    END IF
  END DO

!-------------------------------------------------------------------------------------------------!
!  2. modification : fixed integral depth to depth array
!  nmodel_unit+ depth_seg:  from subrountine load_integral_depth
!  and ismulti is alway F for nmodel_unit is fixed

!%------------------------------------------------------------%
  !| 2.# of model grids to be inverted (equally spacing)          |
  !| For multi-scale parameterization, # of nmodel_unit has to  |
  !| be a poer of 2 plus 1                                      |
  !%------------------------------------------------------------%

  WRITE(*,'(A,I7,A,F10.1)') ' # of save grids=', nmodel_unit, &
                            ', truncate depth (km)=', truc_depm/km2m
  !%------------------------------------------------%
  !| Find the index of radius, such that            |
  !|   model_1d(i1)%rr <= rbound(2)                 |
  !|   model_1d(i2)%rr >= rbound(1)                 |
  !|   model_1d(i1)%rr < model_1d(i2)%rr        |
  !%------------------------------------------------%
  ! NOTE:
  !----------------------------------------------------------------------
! Given an array xx(1:N), and given a value x, returns a value j such
! that x is between xx(j) and xx(j+1). xx must be monotonic, either
! increasing or decreasing. j = 0 or j = N is returned to indicate
! that x is out of range.
! mean: when x is radius of earth, sub locate will return N-1 rather
! than N
!---------------------------------------------------------------------

  IF (ALLOCATED(irbds)) DEALLOCATE(irbds)
  ALLOCATE(irbds(3,nmodel_unit),STAT=Ier)
  IF (Ier/=0) STOP 'allocate error of irbds in open_fmerge_gram'

  DO i = 1, nmodel_unit !\ to greater depth

    rbound_up=rn-depth_seg(i)%up
    rbound_dw=rn-depth_seg(i)%dw
    irbds(1,i)=locate(model_1d%rr,rbound_dw)
    irbds(2,i)=locate(model_1d%rr,rbound_up)
    cdep=(rbound_up+rbound_dw)/2_r8b
    irbds(3,i)=locate(model_1d%rr,cdep)

     IF (rbound_up==model_1d(n)%rr) irbds(2,i)=n  !/locate mis

  END DO



!-------------------------------------------------------------------------------------------------- !
  !%-----------------------------------------------%
  !| if multi-domain parameterization is requiared |
  !%-----------------------------------------------%
!  nmodel=0
!  IF ( isbeta  ) nmodel=nmodel+nmodel_unit
!  IF ( isalpha ) nmodel=nmodel+nmodel_unit
!  IF ( isrho   ) nmodel=nmodel+nmodel_unit
!
!  !%-------------------------------------------------------%
!  !| allocate gram_row and open fmerge      |
!  !%-------------------------------------------------------%
!  IF ( ALLOCATED(gram_row) ) DEALLOCATE(gram_row)
!  ALLOCATE(gram_row(nmodel),STAT=Ier)
!  IF (Ier/=0) STOP 'allocate error of gram_row in sub open_fmerge_gram'
!  INQUIRE(IOLENGTH=iolenGde) gram_row, data0(1)%dd, data0(1)%derr
!
!
!  OPEN(UNIT=ifid_merge, FILE=fmerge,FORM='UNFORMATTED',ACCESS='DIRECT',&
!       ACTION='WRITE',RECL=iolenGde,IOSTAT=Ier)
!   IF (Ier/=0) THEN
!     CLOSE(UNIT=ifid_merge,STATUS='DELETE')
! !    WRITE(*,*) 'Close file:', TRIM(fmerge)
!     OPEN(UNIT=ifid_merge, FILE=fmerge,FORM='UNFORMATTED',ACCESS='DIRECT',&
!          ACTION='WRITE',RECL=iolenGde,IOSTAT=Ier)
!   END IF


  !%----------------------------------%
  !| write out model information file |
  !%----------------------------------%
!  fmergepp=TRIM(fmerge)//'pp'
!  INQUIRE(FILE=fmergepp,EXIST=isexist)
!  IF ( .NOT. isexist ) THEN
!    OPEN(UNIT=fid1,FILE=fmergepp,ACTION='WRITE')
!      WRITE(UNIT=fid1, FMT='(3I12,E20.7,E20.8)') ndata, nmodel, nmodel_unit, &
!           minfo%xb, minfo%dx
!      WRITE(UNIT=fid1,FMT='(A/A)') 'T', 'T'  !\ ismks, isderr
!      WRITE(UNIT=fid1,FMT=*) is_relative_perturb
!      WRITE(UNIT=fid1,FMT=*) isbeta, isalpha, isrho
!      IF ( isbeta ) &
!        WRITE(UNIT=fid1,FMT='(E20.7)') (model_1d(irbds(3,i))%vs,i=1,nmodel_unit)
!      IF ( isalpha ) &
!        WRITE(UNIT=fid1,FMT='(E20.7)') (model_1d(irbds(3,i))%vp,i=1,nmodel_unit)
!      IF ( isrho ) &
!        WRITE(UNIT=fid1,FMT='(E20.7)') (model_1d(irbds(3,i))%rho,i=1,nmodel_unit)
!      WRITE(UNIT=fid1,FMT='(2E16.7)') (model_1d(irbds(3,i))%Qa, model_1d(irbds(3,i))%Qmu,i=1,nmodel_unit)
!    CLOSE(UNIT=fid1)
!  END IF

  !%-----------------------------------------------%
  !| 3. change the  direct write gram_row to save using dcdC,dcdA,dcdL
  !%-----------------------------------------------%


!  IF ( ALLOCATED(dcdC) .OR. ALLOCATED(dcdA) .OR. ALLOCATED(dcdL) ) STOP 'wrong only S mode 0 branch'
!! ndata: all input data number
!  ALLOCATE(dcdC(nmodel_unit,ndata),dcdA(nmodel_unit,ndata),dcdL(nmodel_unit,ndata),STAT=Ier)
!  IF(Ier/=0) STOP'allocate error of dcdC,dcdA,dcdL in sub open_fmerge_gram'



  END SUBROUTINE open_fmerge_gram

  SUBROUTINE merge_kern(ffre,fkern,fmergeid,ftxtid,f_target,&
                        isphs,idata,isfalse_kern,NNL,dcdL,dcdA,dcdC,dcdF,nmodel_unit)
  ! NOTE: deleted is_relative_perturb in parameter list
  ! -> dcdL,dcdA
  ! use isalpha, isbeta, isrho, nhead, irbds, model_1d,&
  !     rn, rhon, wn2
  !---------------------------------------------------------------
  ! read mineos_kernel binary output file and then
  ! merge kernel, data misfit (also data error) in to a binaray
  ! from a trucation depth to surface
  !---------------------------------------------------------------
  CHARACTER(LEN=*), INTENT(IN) :: ffre, fkern
  INTEGER, INTENT(IN) :: fmergeid, ftxtid, idata
  REAL(KIND=r8b), INTENT(IN) :: f_target !\in mhz
  LOGICAL, INTENT(IN) :: isphs !\is data phase velocity
!  LOGICAL, INTENT(IN) :: is_relative_perturb !\ is relative or absolute perturbation kernel
  LOGICAL, INTENT(OUT) :: isfalse_kern  !\ succesful calculate kernel or not

  INTEGER, PARAMETER :: fid1=31
  INTEGER :: imode, nr_unit
  REAL(KIND=r8b), ALLOCATABLE :: kern( : )
  REAL(KIND=r8b) :: pdata
  LOGICAL :: isexist
  INTEGER :: Ier, ishift, nsave, ib, ie, ic, i
  INTEGER :: mrec
 !% Kernel part
  REAL(KIND=r8b), ALLOCATABLE :: kernC( : ),kernA(:),kernL(:),kernF(:)
!  REAL(KIND=r8b), ALLOCATABLE :: kernN(:),kernF(:),kernRho(:),kernD(:)
 LOGICAL :: isCAL=.TRUE.
!% Norm part
   real*8 dnorm,tnorm,fnorm,gravity,vnorm,Lnorm,rhonorm
   REAL(KIND=r8b) ::dcdvn, dcdpn, dcddn, dcdetan, rnn ,dcdLn
   data dnorm,gravity/5515.d0,6.6723d-11/
   REAL*8 zero, eps, eps2
   REAl*8 pi
! Added by C. M. Liu
   INTEGER,INTENT(IN):: NNL
   REAL(KIND=r8b) dcdA(NNL,NNL),dcdL(NNL,NNL)
   REAL(KIND=r8b) dcdC(NNL,NNL),dcdF(NNL,NNL)
   INTEGER,INTENT(IN) :: nmodel_unit
!   rn=6371000.d0
  zero=0.d0
  eps =1.d-37
  eps2=EPSILON(1.d0)*1E+4
  pi=4.d0*datan(1.d0)
  tnorm=1.d0/dsqrt(pi*dnorm*gravity)
  fnorm=1.d0/tnorm

  isfalse_kern = .FALSE.
  INQUIRE(FILE=fkern,EXIST=isexist)
  IF ( .NOT. isexist ) THEN
    WRITE(*,*) 'kernel: ', idata, 'does not exist. Return!'
    isfalse_kern = .TRUE.
    RETURN
  END IF

  !%----------------------------------------------%
  !| grep the proper mode from mineos ASCII file  |
  !%----------------------------------------------%
  imode=get_imode(ffre,f_target)

  !%-----------------------------------------------------------%
  !| read kernel from binary file and store in variable: kern  |
  !|  -> kern (total length: nhead+4*nr_unit)                  |
  !%-----------------------------------------------------------%
  CALL read_mineos_kern(fkern,imode,isphs,pdata,nr_unit,mrec)
  !/ idata, n, l, fhz, Qmod, err, pdata(pv or gv) , data, derr
  !   (vel is in km/s)
  WRITE(UNIT=ftxtid,FMT='(I5,1X,I3,3(1X,E14.7),1X,E14.6,3(1X,F9.4))')&
    idata, NINT(kern(1)), kern(2), kern(3)/twopi, kern(4), kern(6), &
    pdata/km2m, data0(idata)%dd/km2m, data0(idata)%derr/km2m

  !%--------------------------------------------------------------%
  !| write useful kernel into gram_row.                           |
  !| Note that gram_row=kern*grid volume                          |
  !| For dc=int(Km dm,dr), Km is kernel                           |
  !| where m can be alpha, beta, density, Therefore unit of these |
  !| kernels (if mks in used)                                     |
  !|   K_alpha: 1/m                                               |
  !|   K_beta : 1/m                                               |
  !|   K_density: m^3/g-s                                         |
  !| if length is in m and kern (mineos output) is in             |
  !| non-dimensional, kernels should be scaled by normalization   |
  !| factor:                                                      |
  !|   Kern_v/rn                                                  |
  !|   Kern_density*wn/rhon                                       |
  !| if length is in km and weight in kg, then                    |
  !|   Kern_v/rn*1000                                             |
  !|   Kern_density*wn/rhon*1000                                  |
  !| Since velocity and density kernels in mks unit have identical|
  !| dynamic range while those in km-kgs unit have over 6-order   |
  !| deviation, I select mks unit                                 |
  !|                                                              |
  !| Replacing dm, Km with dm/m and m*Km respectively,            |
  !| m*Km is regarded as "relative perturbation kernel"           |
  !%--------------------------------------------------------------%

!% -----------------------------------------------------------------------------------------------------------%!
!%   Normalization
!%   Kernel order: A,C,L,N.F.rho,d
!%------------------------------------------------------------------------------------------------------------%!
!%  m-kg-s unite
   dcdvn = 1.d0/rn         ! 1/m
   dcdpn=fnorm/dnorm       ! m^3/kg*1/s
   rnn=1.d0
   dcddn=fnorm
   dcdLn=tnorm/(dnorm*rn**2) ! ms/kg
   dcdetan=fnorm
!%convert m-kg-s to km-g/cm^3-s
   dcdLn=dcdLn*1.d9
   dcdpn=dcdpn*1.d3
   dcdvn=dcdvn*1.d3
!% -----------------------------------------------------------------------------------------------------------%!
!%   Kernel output
!%   Kernel order: C,A,L,N.F.rho,d
!%------------------------------------------------------------------------------------------------------------%!
    ALLOCATE(kernC(nr_unit),kernA(nr_unit),kernL(nr_unit) ,kernF(nr_unit),STAT=Ier)
    IF (Ier/=0) STOP 'allocate error in sub merge_kern, Terminted!'

!    ALLOCATE(kernN(nr_unit),kernF(nr_unit) ,STAT=Ier)
!    IF (Ier/=0) STOP 'allocate error in sub merge_kern, Terminted!'

 do i=nhead+1,mrec
          kern(i-nhead)=kern(i)
  end do
  DO i=1,nr_unit
       kernC(i)=dble(kern(i))*dcdLn        !dc/dC
       kernA(i)=dble(kern(nr_unit+i))*dcdLn      !dc/dA
       kernL(i)=dble(kern(2*nr_unit+i))*dcdLn    !dc/dL
       kernF(i)=dble(kern(4*nr_unit+i))*dcdLn    !dc/dF
!       kernN(i)=dble(kern(3*nr_unit+i))*dcdLn    !dc/dN
!       kernRho(i)=dble(kern(5*nr_unit+i))*dcdpn    !dc/drho
!       kernD(i)=dble(kern(6*nr_unit+i))*dcddn    !dc/dd

!   Setting numerically nonzero but insignificantly small functional values
!   to exactly zero to avoid later numerical problems.
       if(dabs(kernC(i)).le.eps) kernC(i)=zero
       if(dabs(kernA(i)).le.eps) kernA(i)=zero
       if(dabs(kernL(i)).le.eps) kernL(i)=zero
        if(dabs(kernF(i)).le.eps) kernF(i)=zero
!       if(dabs(kernN(i)).le.eps) kernN(i)=zero
!      if(dabs(kernRho(i)).le.eps) kernRho(i)=zero
!      if(dabs(kernD(i)).le.eps) kernD(i)=zero
  ENDDO


!%  C=p*Vpv^2
!% A=p*Vph^2
!% L=p*Vsv^2
    IF ( isCAL ) THEN
!  irbds: from surface to deeper, depth ascending order.
!  kernA, model_1d: in radius ascending order, depth descending order.
      DO i = 1, nmodel_unit
        ib=irbds(1,i)
        ie=irbds(2,i)
        dcdC(idata,i)=integ(model_1d(ib:ie)%rr, kernC(ib:ie))*m2km
        dcdA(idata,i)=integ(model_1d(ib:ie)%rr, kernA(ib:ie))*m2km
        dcdL(idata,i)=integ(model_1d(ib:ie)%rr, kernL(ib:ie))*m2km
        dcdF(idata,i)=integ(model_1d(ib:ie)%rr,kernF(ib:ie))*m2km
      END DO
  END IF

  DEALLOCATE(kern)
  DEALLOCATE(kernC,kernA,kernL,kernF)


  CONTAINS
    FUNCTION get_imode(ffre,f_target)
    !%------------------------------------------------------------%
    !| find the ith of mode such that its frequency is closest to |
    !| observation.                                               |
    !%------------------------------------------------------------%
    CHARACTER(LEN=*), INTENT(IN) :: ffre  !\mineos output ASCII file
    REAL(KIND=r8b),INTENT(IN) :: f_target !\target mode frequency (hz)
    INTEGER :: get_imode

    INTEGER, PARAMETER :: fid2=51
    REAL(KIND=r8b) :: ftbl, dfmhz, dfmhzmin, r8btmp
    INTEGER :: i, Ier, itmp
    CHARACTER(LEN=2) :: cjcom

    get_imode=0
    OPEN(UNIT=fid2,FILE=ffre,STATUS='OLD',ACTION='READ')
    dfmhzmin=HUGE(1.0_r8b)
    i=0
    DO
      READ(UNIT=fid2,FMT=*,IOSTAT=Ier) &
      itmp, cjcom, r8btmp, r8btmp, ftbl, r8btmp, r8btmp, r8btmp, r8btmp
      IF (Ier<0) EXIT
      i=i+1
      dfmhz=ABS(ftbl-f_target)
      IF (dfmhz<dfmhzmin) THEN
        get_imode=i
        dfmhzmin=dfmhz
      END IF
    END DO
    CLOSE(fid2)
    END FUNCTION get_imode

    SUBROUTINE read_mineos_kern(fkern,imode,isphs,pdata,nr_unit,mrec)
    !-> kern, model_1d
    ! USE ModeNorm_module, ONLY: wn2, rn
    ! use nrmax, nhead, km2m
    !--------------------------------------------------------
    ! read kernel from binary file and store in variable kern
    !--------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: fkern !\ mineos output kernel file
    INTEGER, INTENT(IN) :: imode !\ the ith mode
    LOGICAL, INTENT(IN) :: isphs !\ is data phase velocity
    REAL(KIND=r8b), INTENT(OUT) :: pdata !\ predicted data
    INTEGER, INTENT(OUT) :: nr_unit !\ amount of saved kernel grids

    INTEGER, PARAMETER :: fid2=50

    REAL(KIND=r8b) :: pv, gv
    INTEGER :: mrec, nr, jcom, n, nic, noc, nhead1
    INTEGER :: Ier, j, i, itmp
    REAL    :: rtmp, tdum(nrmax,9)
    LOGICAL :: ismatch

    nhead1=nhead+1
    OPEN(UNIT=fid2,FILE=fkern,STATUS='OLD',ACTION='READ',&
         FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    !%------------------------------------------------------------%
    !| Read the header of kernel file for info on reference model |
    !| for Toroidal mode: saved knots is from noc+1 to surface    |
    !%------------------------------------------------------------%
!    READ(UNIT=fid2) jcom, wmin, wmax, flmin, flmax, wgrav
!    READ(UNIT=fid2) n,nic,noc,ifanis,trefl,((tdum(i,j),i=1,n),j=1,9)
    READ(UNIT=fid2) jcom, rtmp, rtmp, rtmp, rtmp, rtmp
    READ(UNIT=fid2) n, nic, noc, itmp, rtmp, ((tdum(i,j),i=1,n),j=1,9)

    !%-----------------------%
    !| get nr, mrec          |
    !%-----------------------%
 !  nkern=7; nhead=6
    IF(jcom==iflag_sph)THEN
      mrec=nkern*n+nhead
      nr=n
    ELSE IF (jcom == iflag_tor) THEN
      mrec=nkern*(n-noc)+nhead
      nr=n-noc
    ELSE
      STOP 'jcom is out of range in sub read_mineos_kern'
    END IF
    ALLOCATE(kern(mrec),STAT=Ier)
    IF (Ier/=0) STOP 'allocate error in sub read_mineos_kern, Terminted!'

    kern=0.0_r8b
    ismatch=.FALSE.
    i=1
    DO
      READ(UNIT=fid2, IOSTAT=Ier) (kern(j),j=1,mrec)
      ismatch = (i==imode)
      IF (ismatch) THEN
       !%--------------------------------%
       !| find the match mode            |
       !%--------------------------------%
        !/ phase velocity in m/s
        pv=kern(3)/SQRT(kern(2)*(kern(2)+1.0_r8b))*rn  !kern(3):angle frequence(2pi*f, f: Hz)   kern(2):l angle order
        !/ group velocity in m/s
        gv=kern(5)*km2m
!        WRITE(*,*) pv**2/gv/kern(3)/km2m
      END IF
      IF (Ier<0.OR.ismatch) EXIT
      i=i+1
    END DO
    CLOSE(UNIT=fid2)

    IF ( isphs ) THEN
      pdata=pv
    ELSE
      pdata=gv
    END IF
    nr_unit=nr

    END SUBROUTINE read_mineos_kern
  END SUBROUTINE merge_kern



  FUNCTION integ(x,y)

  !--------------------------------------------------------------------
  ! get integral of y(x) by trapzoidal method
  !-------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8b) :: integ
    REAL(KIND=r8b), DIMENSION( : ), INTENT(IN) :: x, y

    REAL(KIND=r8b) :: dx
    INTEGER :: n, n1, i

    n=SIZE(x)
    n1=n-1
    IF (n/=SIZE(y)) STOP 'size of x and y has to be the same in sub int...'
    integ =0.0_r8b
    DO i = 1, n1
      dx=x(i+1)-x(i)
      integ=integ+SUM(y(i:i+1))*dx
    END DO
    integ=integ/2.0_r8b
  END FUNCTION integ

END MODULE buildG_dcdm1D_module
