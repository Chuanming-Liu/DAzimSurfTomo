MODULE sort_module
!----------------------------------------------------------------------
! Adaption of "Numerical recipes" (fortran90, v210)
!---------------------------------------------------------------------
  USE data_type
  USE nrutil, ONLY : arth,assert_eq,nrerror,swap 

  INTERFACE indexx
    MODULE PROCEDURE indexx_sp
    MODULE PROCEDURE indexx_i4b
  END INTERFACE indexx

  PRIVATE
  PUBLIC :: indexx

CONTAINS
	SUBROUTINE indexx_sp(arr,indexa)
	IMPLICIT NONE
	REAL(KIND=r4b), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(KIND=i4b), DIMENSION(:), INTENT(OUT) :: indexa
	INTEGER(KIND=i4b), PARAMETER :: NN=15, NSTACK=50
	REAL(KIND=r4b) :: a
	INTEGER(KIND=i4b) :: n,k,i,j,indext,jstack,l,r
	INTEGER(KIND=i4b), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(indexa),size(arr),'indexx_sp')
	indexa=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=indexa(j)
				a=arr(indext)
				do i=j-1,l,-1
					if (arr(indexa(i)) <= a) exit
					indexa(i+1)=indexa(i)
				end do
				indexa(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(indexa(k),indexa(l+1))
			call icomp_xchg(indexa(l),indexa(r))
			call icomp_xchg(indexa(l+1),indexa(r))
			call icomp_xchg(indexa(l),indexa(l+1))
			i=l+1
			j=r
			indext=indexa(l+1)
			a=arr(indext)
			do
				do
					i=i+1
					if (arr(indexa(i)) >= a) exit
				end do
				do
					j=j-1
					if (arr(indexa(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(indexa(i),indexa(j))
			end do
			indexa(l+1)=indexa(j)
			indexa(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(KIND=i4b), INTENT(INOUT) :: i,j
	INTEGER(KIND=i4b) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_sp

	SUBROUTINE indexx_i4b(iarr,indexa)
	IMPLICIT NONE
	INTEGER(KIND=i4b), DIMENSION(:), INTENT(IN) :: iarr
	INTEGER(KIND=i4b), DIMENSION(:), INTENT(OUT) :: indexa
	INTEGER(KIND=i4b), PARAMETER :: NN=15, NSTACK=50
	INTEGER(KIND=i4b) :: a
	INTEGER(KIND=i4b) :: n,k,i,j,indext,jstack,l,r
	INTEGER(KIND=i4b), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(indexa),size(iarr),'indexx_i4b')
	indexa=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=indexa(j)
				a=iarr(indext)
				do i=j-1,l,-1
					if (iarr(indexa(i)) <= a) exit
					indexa(i+1)=indexa(i)
				end do
				indexa(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(indexa(k),indexa(l+1))
			call icomp_xchg(indexa(l),indexa(r))
			call icomp_xchg(indexa(l+1),indexa(r))
			call icomp_xchg(indexa(l),indexa(l+1))
			i=l+1
			j=r
			indext=indexa(l+1)
			a=iarr(indext)
			do
				do
					i=i+1
					if (iarr(indexa(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(indexa(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(indexa(i),indexa(j))
			end do
			indexa(l+1)=indexa(j)
			indexa(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(KIND=i4b), INTENT(INOUT) :: i,j
	INTEGER(KIND=i4b) :: swp
	if (iarr(j) < iarr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_i4b

END MODULE sort_module
