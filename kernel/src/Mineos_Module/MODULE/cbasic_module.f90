MODULE data_type
  IMPLICIT NONE
  INTEGER, PARAMETER :: r8b = SELECTED_REAL_KIND( p=15, r=300 ),&
                        i8b = SELECTED_INT_KIND ( r=18 ),&
                        r4b = SELECTED_REAL_KIND( p=6, r=30 ), &
                        i4b = SELECTED_INT_KIND ( r=9 ), &
                        i2b = SELECTED_INT_KIND ( r=4 ), &
                        i1b = SELECTED_INT_KIND ( r=2 ), &
                        LGT = KIND(.TRUE.),&
                        r8bc =SELECTED_REAL_KIND( p=15 )  
END MODULE data_type

MODULE math_par_8b
  USE data_type
  IMPLICIT NONE
  REAL(KIND=r8b), PARAMETER :: pi_r8 = 3.141592653589793_r8b, &
                               deg2rad_r8 = 0.01745329251994330_r8b, &
                               rad2deg_r8 = 57.2957795130823_r8b

END MODULE math_par_8b

MODULE math_par_4b
  USE data_type
  IMPLICIT NONE
  REAL, PARAMETER :: pi = 3.1415926, &
                     deg2rad = 0.01745329252, &
                     rad2deg = 57.29578

END MODULE math_par_4b


