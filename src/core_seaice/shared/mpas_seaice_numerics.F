!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_numerics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_numerics

  use mpas_derived_types
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_solve_linear_basis_system

contains

!-----------------------------------------------------------------------
! LU decomposition
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_solve_linear_basis_system
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         leftMatrix
    real(kind=RKIND), dimension(:), intent(in) :: &
         rightHandSide
    real(kind=RKIND), dimension(:), intent(out) :: &
         solutionVector

    real(kind=RKIND), allocatable, dimension(:,:) :: a
    real(kind=RKIND), allocatable, dimension(:) :: b
    integer, allocatable, dimension(:) :: indices
    real(kind=RKIND) :: d

    integer :: n

    n = size(rightHandSide)

    allocate(a(n,n))
    allocate(b(n))
    allocate(indices(n))

    a(:,:) = leftMatrix(:,:)
    b(:) = rightHandSide(:)

    call lu_decomposition(a,indices,d)

    call lu_back_subsitution(a,indices,b)

    solutionVector(:) = b(:)

    deallocate(a)
    deallocate(b)
    deallocate(indices)

  end subroutine seaice_solve_linear_basis_system

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  lu_decomposition
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-------------------------------------------------------------

  subroutine lu_decomposition(a,indices,d)

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         a

    integer, dimension(:), intent(out) :: &
         indices

    real(kind=RKIND), intent(out) :: &
         d

    real(kind=RKIND), dimension(size(a,1)) :: &
         maxa

    real(kind=RKIND), parameter :: &
         tiny = 1.0e-20_RKIND

    integer :: &
         j, &
         n, &
         jmax

    ! check sizes
    if (size(a,1) == size(a,2) .and. size(a,2) == size(indices)) then
       n = size(a,1)
    else
       call mpas_log_write("lu_decomposition: unequal array sizes: $i $i $i", &
            MPAS_LOG_CRIT, intArgs=(/size(a,1),size(a,2),size(indices)/))
    endif

    d = 1.0_RKIND
    maxa = maxval(abs(a),dim=2)
    if (any(maxa == 0.0_RKIND)) then
       call mpas_log_write("lu_decomposition: singular matrix", &
            MPAS_LOG_CRIT)
    endif
    maxa = 1.0_RKIND / maxa

    do j = 1, n

       jmax = (j-1) + get_array_max_location(maxa(j:n)*abs(a(j:n,j)))

       if (j /= jmax) then
          call swap_arrays(a(jmax,:),a(j,:))
          d = -d
          maxa(jmax) = maxa(j)
       end if

       indices(j) = jmax
       if (a(j,j) == 0.0_RKIND) &
            a(j,j) = tiny
       a(j+1:n,j)     = a(j+1:n,j)     / a(j,j)
       a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - outer_product(a(j+1:n,j),a(j,j+1:n))

    end do

  end subroutine lu_decomposition

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  lu_back_subsitution
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-------------------------------------------------------------

  subroutine lu_back_subsitution(a,indices,b)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         a

    integer, dimension(:), intent(in) :: &
         indices

    real(kind=RKIND), dimension(:), intent(inout) :: &
         b

    integer :: &
         i, &
         n, &
         j, &
         k

    real(kind=RKIND) :: &
         sums

    if (size(a,1) == size(a,2) .and. size(a,2) == size(indices)) then
       n = size(a,1)
    else
       call mpas_log_write("lu_back_subsitution: unequal array sizes: $i $i $i", &
            MPAS_LOG_CRIT, intArgs=(/size(a,1),size(a,2),size(indices)/))
    endif

    j = 0

    do i = 1, n

       k = indices(i)
       sums = b(k)
       b(k) = b(i)

       if (j /= 0) then
          sums = sums - dot_product(a(i,j:i-1), b(j:i-1))
       else if (sums /= 0.0_RKIND) then
          j = i
       end if

       b(i) = sums

    end do

    do i = n, 1, -1
       b(i) = (b(i) - dot_product(a(i,i+1:n), b(i+1:n))) / a(i,i)
    end do

  end subroutine lu_back_subsitution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_array_max_location
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-------------------------------------------------------------

  function get_array_max_location(array) result(max_loc)

    real(kind=RKIND), dimension(:), intent(in) :: &
         array

    integer :: &
         max_loc

    integer, dimension(1) :: &
         max_loc_array

    max_loc_array = maxloc(array(:))
    max_loc       = max_loc_array(1)

  end function get_array_max_location

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  outer_product
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-------------------------------------------------------------

  function outer_product(x, y) result(outerproduct)

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, &
         y

    real(kind=RKIND), dimension(size(x),size(y)) :: &
         outerproduct

    outerproduct = &
         spread(x,dim=2,ncopies=size(y)) * &
         spread(y,dim=1,ncopies=size(x))

  end function outer_product

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  swap_arrays
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-------------------------------------------------------------

  subroutine swap_arrays(x, y)

    real(kind=RKIND), dimension(:), intent(inout) :: &
         x, &
         y

    real(kind=RKIND), dimension(size(x)) :: &
         tmp

    tmp = x
    x   = y
    y   = tmp

  end subroutine swap_arrays

  !-------------------------------------------------------------

end module seaice_numerics
