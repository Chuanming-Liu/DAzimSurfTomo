MODULE interpol_module
  USE nrutil, ONLY : nrerror, assert_eq
  USE data_type

  INTERFACE splint
    MODULE PROCEDURE splint_scalar
    MODULE PROCEDURE splint_ary
  END INTERFACE splint

  PRIVATE
  PUBLIC :: spline, splint, locate, tridag_par, polint

CONTAINS

  SUBROUTINE spline(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: x,y
    REAL(KIND=r8b), INTENT(IN) :: yp1,ypn
    REAL(KIND=r8b), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(KIND=i4b) :: n
    REAL(KIND=r8b), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_r8b*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_r8b*(c(2:n-1)+a(2:n-1))
    b(1)=1.0_r8b
    b(n)=1.0_r8b
    if (yp1 > 0.99E300_r8b) then
        r(1)=0.0_r8b
        c(1)=0.0_r8b
    else
        r(1)=(3.0_r8b/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        c(1)=0.5_r8b
    end if
    if (ypn > 0.99E300_r8b) then
        r(n)=0.0_r8b
        a(n)=0.0_r8b
    else
        r(n)=(-3.0_r8b/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
        a(n)=0.5_r8b
    end if
    call tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline

  FUNCTION splint_scalar(xa,ya,y2a,x)
!----------------------------------------------------------------------
! Given the arrays xa and ya, which tabulates a function (with the
! xa_i's in increasing or decreasing order), and given the array
! y2a, which is the output from another SUBROUTINE "SPLINE" and
! represents the second derivative of ya (i.e. y2a = d2ya/dxa2), and
! given a value of x, and then this routine returns a cubic-spline
! interpolated value. The arrays xa, ya, y2a are all of the same
! size.
!---------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(KIND=r8b), INTENT(IN) :: x
    REAL(KIND=r8b) :: splint_scalar
    INTEGER(KIND=i4b) :: khi,klo,n
    REAL(KIND=r8b) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0_r8b) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_scalar=a*ya(klo)+b*ya(khi)+&
                  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_r8b
  END FUNCTION splint_scalar

  FUNCTION splint_ary(xa,ya,y2a,x)
!----------------------------------------------------------------------
! similar to splint, but intput "x" is a 1D array
!----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: x
    REAL(KIND=r8b), DIMENSION(SIZE(x)) :: splint_ary
    INTEGER(KIND=i4b) :: n, nx, i
    INTEGER(KIND=i4b), DIMENSION(SIZE(x)) :: khi,klo
    REAL(KIND=r8b), DIMENSION(SIZE(x)) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    nx=SIZE(x)
    DO i = 1, nx
       klo(i)=max(min(locate(xa,x(i)),n-1),1)
    END DO
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (ANY(h == 0.0_r8b)) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_ary=a*ya(klo)+b*ya(khi)+&
               ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_r8b
  END FUNCTION splint_ary

  FUNCTION locate(xx,x)
!----------------------------------------------------------------------
! Given an array xx(1:N), and given a value x, returns a value j such
! that x is between xx(j) and xx(j+1). xx must be monotonic, either
! increasing or decreasing. j = 0 or j = N is returned to indicate
! that x is out of range.
!---------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: xx
    REAL(KIND=r8b), INTENT(IN) :: x
    INTEGER(KIND=i4b) :: locate
    INTEGER(KIND=i4b) :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
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

  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    IMPLICIT NONE
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(KIND=r8b), DIMENSION(:), INTENT(OUT) :: u
    INTEGER(KIND=i4b), PARAMETER :: NPAR_TRIDAG=4
    INTEGER(KIND=i4b) :: n,n2,nm,nx
    REAL(KIND=r8b), DIMENSION(size(b)/2) :: y,q,piva
    REAL(KIND=r8b), DIMENSION(size(b)/2-1) :: x,z
    REAL(KIND=r8b), DIMENSION(size(a)/2) :: pivc
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
    if (n < NPAR_TRIDAG) then
            call tridag_ser(a,b,c,r,u)
    else
            if (maxval(abs(b(1:n))) == 0.0_r8b) &
                    call nrerror('tridag_par: possible singular matrix')
            n2=size(y)
            nm=size(pivc)
            nx=size(x)
            piva = a(1:n-1:2)/b(1:n-1:2)
            pivc = c(2:n-1:2)/b(3:n:2)
            y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
            q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
            if (nm < n2) then
                    y(n2) = b(n)-piva(n2)*c(n-1)
                    q(n2) = r(n)-piva(n2)*r(n-1)
            end if
            x = -piva(2:n2)*a(2:n-2:2)
            z = -pivc(1:nx)*c(3:n-1:2)
            call tridag_par(x,y,z,q,u(2:n:2))
            u(1) = (r(1)-c(1)*u(2))/b(1)
            u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
                    -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
            if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
    CONTAINS
      SUBROUTINE tridag_ser(a,b,c,r,u)
        IMPLICIT NONE
        REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: a,b,c,r
        REAL(KIND=r8b), DIMENSION(:), INTENT(OUT) :: u
        REAL(KIND=r8b), DIMENSION(size(b)) :: gam
        INTEGER(KIND=i4b) :: n,j
        REAL(KIND=r8b) :: bet
        n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
        bet=b(1)
        if (bet == 0.0_r8b) call nrerror('tridag_ser: Error at code stage 1')
        u(1)=r(1)/bet
        do j=2,n
                gam(j)=c(j-1)/bet
                bet=b(j)-a(j-1)*gam(j)
                if (bet == 0.0_r8b) &
                        call nrerror('tridag_ser: Error at code stage 2')
                u(j)=(r(j)-a(j-1)*u(j-1))/bet
        end do
        do j=n-1,1,-1
                u(j)=u(j)-gam(j+1)*u(j+1)
        end do
        END SUBROUTINE tridag_ser
  END SUBROUTINE tridag_par

  SUBROUTINE polint(xa,ya,x,y,dy)
    REAL(KIND=r8b), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(KIND=r8b), INTENT(IN) :: x
    REAL(KIND=r8b), INTENT(OUT) :: y,dy
    INTEGER(KIND=i4b) :: m,n,ns
    REAL(KIND=r8b), DIMENSION(size(xa)) :: c,d,den,ho
    !if size(xa)==size(yz),assign size to n
    n=assert_eq(size(xa),size(ya),'polint') !if size(xa)==size(yz)
    c=ya
    d=ya
    ho=xa-x
    ns=MINLOC(ABS(x-xa),DIM=1)
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
            den(1:n-m)=ho(1:n-m)-ho(1+m:n)
            if (any(den(1:n-m) == 0.0)) then
                write(*,*)'xa array (f mHz):', (xa(i),i=1,n)
                write(*,*)'ya array (angular l):', (ya(i),i=1,n)
                write(*,*)'target f Hz',x
                call nrerror('polint: calculation failure')
           end if
            den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
            d(1:n-m)=ho(1+m:n)*den(1:n-m)
            c(1:n-m)=ho(1:n-m)*den(1:n-m)
            if (2*ns < n-m) then
                    dy=c(ns+1)
            else
                    dy=d(ns)
                    ns=ns-1
            end if
            y=y+dy
    end do
  END SUBROUTINE polint

END MODULE interpol_module
