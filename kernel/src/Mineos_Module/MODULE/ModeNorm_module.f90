MODULE modenorm_module
  USE data_type
  USE math_par_8b, ONLY : pi_r8

  IMPLICIT NONE
  REAL(KIND=r8b), PARAMETER :: &
                     bigg=6.6723e-11_r8b, &
                     rhon=5515.0_r8b, &
                     rn=6371000.0_r8b, &
                     wn2=pi_r8*bigg*rhon, &
                     gn =wn2*rn
!                     wn=SQRT(wn2), &
!                     vn=rn*wn

  PRIVATE
!  PUBLIC :: rhon, rn, gn, wn, vn, bigg, wn2
  PUBLIC :: rhon, rn, gn, bigg, wn2
END MODULE modenorm_module
