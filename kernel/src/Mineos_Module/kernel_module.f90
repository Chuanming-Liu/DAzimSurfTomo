MODULE kernel_module
!    NLAY  - layers in model: check getmodel.f for the value
!    maxdata: maximun number of observed data
!   dcdA(model,period)
  IMPLICIT NONE
  INTEGER maxdata,NL
  Parameter(maxdata = 200,NL=200)
!  REAL(KIND=8),ALLOCATABLE :: dcdC(:,:),dcdA(:,:),dcdL(:,:)
  REAL(KIND=8) :: dcdA(NL,maxdata), dcdL(NL,maxdata)
END MODULE kernel_module



