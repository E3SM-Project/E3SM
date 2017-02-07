!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_numerics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_numerics

  use mpas_derived_types

  implicit none

  private
  save

  public :: &
       cice_solve_linear_basis_system

contains

!-----------------------------------------------------------------------
! LU decomposition
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_solve_linear_basis_system
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         leftMatrix
    real(kind=RKIND), dimension(:), intent(in) :: &
         rightHandSide
    real(kind=RKIND), dimension(:), intent(out) :: &
         solutionVector

    real(kind=RKIND), allocatable, dimension(:,:) :: a
    real(kind=RKIND), allocatable, dimension(:) :: b
    integer, allocatable, dimension(:) :: indx
    real(kind=RKIND) :: d

    integer :: n

    n = size(rightHandSide)

    allocate(a(n,n))
    allocate(b(n))
    allocate(indx(n))

    a(:,:) = leftMatrix(:,:)
    b(:) = rightHandSide(:)

    call ludcmp(a,indx,d)

    call lubksb(a,indx,b)

    solutionVector(:) = b(:)

    deallocate(a)
    deallocate(b)
    deallocate(indx)

  end subroutine cice_solve_linear_basis_system

  !-------------------------------------------------------------

  SUBROUTINE ludcmp(a,indx,d)

    IMPLICIT NONE
    REAL(kind=RKIND), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT) :: indx
    REAL(kind=RKIND), INTENT(OUT) :: d
    REAL(kind=RKIND), DIMENSION(size(a,1)) :: vv
    REAL(kind=RKIND), PARAMETER :: TINY=1.0e-20_RKIND
    INTEGER :: j,n,imax
    n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d=1.0_RKIND
    vv=maxval(abs(a),dim=2)
    if (any(vv == 0.0_RKIND)) call nrerror('singular matrix in ludcmp')
    vv=1.0_RKIND/vv
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       end if
       indx(j)=imax
       if (a(j,j) == 0.0_RKIND) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    end do
  END SUBROUTINE ludcmp

  !-------------------------------------------------------------

  SUBROUTINE lubksb(a,indx,b)

    IMPLICIT NONE
    REAL(kind=RKIND), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER, DIMENSION(:), INTENT(IN) :: indx
    REAL(kind=RKIND), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n,ii,ll
    REAL(kind=RKIND) :: summ
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0_RKIND) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE lubksb

  !-------------------------------------------------------------

  FUNCTION assert_eq(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq

  !-------------------------------------------------------------

  FUNCTION imaxloc(arr)
    REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc

  !-------------------------------------------------------------

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  !-------------------------------------------------------------

  FUNCTION outerprod(a,b)
    REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: a,b
    REAL(kind=RKIND), DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod

  !-------------------------------------------------------------

  SUBROUTINE swap(a,b)
    REAL(kind=RKIND), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(kind=RKIND), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap

  !-------------------------------------------------------------

end module cice_numerics
