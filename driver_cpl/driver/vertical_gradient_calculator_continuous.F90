module vertical_gradient_calculator_continuous

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type for
  ! computing piecewise continuous vertical gradients using a matrix solve.

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_sys_mod, only : shr_sys_abort
  use shr_matrix_mod, only : tridiagonal_inverse

  implicit none
  private

  public :: vertical_gradient_calculator_continuous_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_continuous_type
     private

     integer :: nelev ! number of elevation classes
     integer :: num_points
     real(r8), allocatable :: field(:,:)   ! field(i,j) is elevation class i, point j
     real(r8), allocatable :: topo(:,:)    ! topo(i,j) is elevation class i, point j

     ! precomputed vertical gradients; vertical_gradient(i,j) is elevation class i, field
     ! j
     real(r8), allocatable :: vertical_gradient(:,:)

     logical, allocatable :: topo_valid(:)  ! whether topo is valid in each point

     ! Bounds of each elevation class. This array has one more element than the number of
     ! elevation classes, since it contains lower and upper bounds for each elevation
     ! class. The indices go 0:nelev. These bounds are guaranteed to be monotonically
     ! increasing.
     real(r8), allocatable :: elevclass_bounds(:)

     ! Calculator to determine initial guesses for gradients. We determine how good the
     ! solution is based on how well we match these initial guesses. We also fall back on
     ! these initial guesses if the gradient in a given elevation class is determined to
     ! be 'bad'.
     class(vertical_gradient_calculator_base_type), allocatable :: calculator_initial_guess

     logical :: calculated  ! whether gradients have been calculated yet

     ! Various statistics for printing diagnostics.
     logical, allocatable :: zeroed_from_topo_out_of_bounds(:) ! [num_points]
     logical, allocatable :: limited_to_zero(:,:) ! [nelev, num_points]
     logical, allocatable :: limited_to_initial_guess(:,:) ! [nelev, num_points]

   contains
     procedure :: calc_gradients
     procedure :: get_gradients_one_class
     procedure :: get_gradients_one_point
     procedure :: print_statistics

     ! This is public so that it can be overridden and/or tested independently by unit
     ! tests
     procedure :: limit_gradients

     procedure, private :: check_topo  ! check topographic heights
     procedure, private :: solve_for_vertical_gradients ! compute vertical gradients for all ECs, for points where we do a matrix solve

     procedure, private :: dl  ! lower half-width of ec
     procedure, private :: du  ! upper half-width of ec
  end type vertical_gradient_calculator_continuous_type

  interface vertical_gradient_calculator_continuous_type
     module procedure constructor
  end interface vertical_gradient_calculator_continuous_type

contains

  !-----------------------------------------------------------------------
  function constructor(field, topo, elevclass_bounds, calculator_initial_guess) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates a vertical_gradient_calculator_continuous_type object.
    !
    ! Pre-condition: elevclass_bounds must be monotonically increasing.
    !
    ! Pre-condition: Topographic heights should all lie inside the bounds of their
    ! respective elevation class (given by elevclass_bounds), with the possible exception
    ! of the lowest elevation class (topographic heights can lie below the arbitrary lower
    ! bound of the elevation class) and the highest elevation class (topographic heights
    ! can lie above the arbitrary upper bound of the elevation class). For grid cells
    ! where this is not true, sets vertical gradient to 0 for all elevation classes.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_continuous_type) :: this  ! function result
    real(r8), intent(in) :: field(:,:)  ! field(i,j) is point i, elevation class j
    real(r8), intent(in) :: topo(:,:)   ! topo(i,j) is point i, elevation class j

    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8)         , intent(in) :: elevclass_bounds(0:)

    ! Initial guesses for gradients. We determine how good the solution is based on how
    ! well we match these initial guesses. We also fall back on these initial guesses if
    ! the gradient in a given elevation class is determined to be 'bad'.
    !
    ! The calc_gradients method should not yet have been called on this object
    class(vertical_gradient_calculator_base_type), intent(in) :: calculator_initial_guess
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%calculated = .false.

    this%num_points = size(field, 1)
    this%nelev = size(field, 2)
    SHR_ASSERT_ALL((ubound(topo) == (/this%num_points, this%nelev/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(elevclass_bounds) == (/this%nelev/)), errMsg(__FILE__, __LINE__))

    allocate(this%elevclass_bounds(0:this%nelev))
    this%elevclass_bounds(:) = elevclass_bounds(:)
    call this%check_elevclass_bounds_monotonic_increasing(this%elevclass_bounds)

    allocate(this%field(this%nelev, this%num_points))
    this%field(:,:) = transpose(field(:,:))
    allocate(this%topo(this%nelev, this%num_points))
    this%topo(:,:) = transpose(topo(:,:))

    allocate(this%topo_valid(this%num_points))
    call this%check_topo()

    allocate(this%vertical_gradient(this%nelev, this%num_points))
    this%vertical_gradient(:,:) = nan

    allocate(this%calculator_initial_guess, source = calculator_initial_guess)

    allocate(this%zeroed_from_topo_out_of_bounds(this%num_points))
    this%zeroed_from_topo_out_of_bounds(:) = .false.
    allocate(this%limited_to_zero(this%nelev, this%num_points))
    this%limited_to_zero(:,:) = .false.
    allocate(this%limited_to_initial_guess(this%nelev, this%num_points))
    this%limited_to_initial_guess(:,:) = .false.

  end function constructor

  !-----------------------------------------------------------------------
  subroutine calc_gradients(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: pt

    character(len=*), parameter :: subname = 'calc_gradients'
    !-----------------------------------------------------------------------

    if (this%calculated) then
       return
    end if

    call this%calculator_initial_guess%calc_gradients()
    do pt = 1, this%num_points
       if (.not. this%topo_valid(pt)) then
          this%vertical_gradient(:,pt) = 0._r8
          this%zeroed_from_topo_out_of_bounds(pt) = .true.
       else
          call this%solve_for_vertical_gradients(pt)
       end if
    end do

    this%calculated = .true.

  end subroutine calc_gradients


  !-----------------------------------------------------------------------
  subroutine get_gradients_one_class(this, elevation_class, gradients)
    !
    ! !DESCRIPTION:
    ! Returns the vertical gradient for all points, at a given elevation class.
    !
    ! this%calc_gradients should already have been called
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_gradients_one_class'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(gradients) == this%num_points), errMsg(__FILE__, __LINE__))

    if (elevation_class < 1 .or. &
         elevation_class > this%nelev) then
       write(logunit,*) subname, ': ERROR: elevation class out of bounds: ', &
            elevation_class, this%nelev
       call shr_sys_abort(subname//': ERROR: elevation class out of bounds')
    end if

    gradients(:) = this%vertical_gradient(elevation_class, :)

  end subroutine get_gradients_one_class

  !-----------------------------------------------------------------------
  subroutine get_gradients_one_point(this, point, gradients)
    !
    ! !DESCRIPTION:
    ! Returns the vertical gradient for all elevation classes, for one point
    !
    ! this%calc_gradients should already have been called
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    integer, intent(in) :: point

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_gradients_one_point'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, errMsg(__FILE__, __LINE__))
    SHR_ASSERT(point <= this%num_points, errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(gradients) == this%nelev), errMsg(__FILE__, __LINE__))

    gradients(:) = this%vertical_gradient(:, point)

  end subroutine get_gradients_one_point

  !-----------------------------------------------------------------------
  subroutine print_statistics(this)
    !
    ! !DESCRIPTION:
    ! Print various statistics on the solve to the logunit
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: ec
    integer :: num_not_out_of_bounds

    character(len=*), parameter :: subname = 'print_statistics'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, errMsg(__FILE__, __LINE__))

    write(logunit, '(a)') "Vertical gradient calculator statistics: "
    write(logunit, '(a, f10.6)') "Fraction with topo out of bounds: ", &
         real(count(this%zeroed_from_topo_out_of_bounds), r8) / real(this%num_points, r8)
    num_not_out_of_bounds = this%num_points - count(this%zeroed_from_topo_out_of_bounds)
    do ec = 1, this%nelev
       write(logunit, '(a, i4, f10.6)') "Remaining fraction limited to 0: ", ec, &
            real(count(this%limited_to_zero(ec,:)), r8) / real(num_not_out_of_bounds, r8)
       write(logunit, '(a, i4, f10.6)') "Remaining fraction limited to initial guess: ", ec, &
            real(count(this%limited_to_initial_guess(ec,:)), r8) / real(num_not_out_of_bounds, r8)
    end do

  end subroutine print_statistics


  !-----------------------------------------------------------------------
  subroutine check_topo(this)
    !
    ! !DESCRIPTION:
    ! Check topographic heights; set this%topo_valid(i) to false if there is a problem in
    ! point i.
    !
    ! Topographic heights in the attribute vector must all lie inside the bounds of their
    ! respective elevation class (given by elevclass_bounds), with the possible exception
    ! of the lowest elevation class (topographic heights can lie below the arbitrary lower
    ! bound of the elevation class) and the highest elevation class (topographic heights
    ! can lie above the arbitrary upper bound of the elevation class)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
    integer :: i

    ! Absolute tolerance for error checks. This is chosen so that it allows for
    ! double-precision roundoff-level errors on values of order 10,000.
    real(r8), parameter :: tol = 1.e-10_r8

    character(len=*), parameter :: subname = 'check_topo'
    !-----------------------------------------------------------------------

    this%topo_valid(:) = .true.

    do i = 1, this%num_points
       do elevclass = 1, this%nelev
          if (elevclass > 1) then
             if (this%topo(elevclass,i) - this%elevclass_bounds(elevclass-1) < -tol) then
                this%topo_valid(i) = .false.
             end if
          end if

          if (elevclass < this%nelev) then
             if (this%topo(elevclass,i) - this%elevclass_bounds(elevclass) > tol) then
                this%topo_valid(i) = .false.
             end if
          end if
       end do
    end do

  end subroutine check_topo

  !-----------------------------------------------------------------------
  subroutine solve_for_vertical_gradients(this, pt)
    !
    ! !DESCRIPTION:
    ! Compute and save vertical gradients for all elevation classes in this point.
    !
    ! This should only be called for points where we have done some initial checks to
    ! show that we should attempt a matrix solve there.
    !
    ! Computes a gradient in each elevation class such that the field is continuous at
    ! interfaces and the sum over squared differences from the mean is minimized.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    integer, intent(in) :: pt  ! point to compute gradients for (1..this%num_points)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: field(this%nelev) ! mean field value of each elevation class
    real(r8) :: topo(this%nelev)  ! mean topo of each elevation class
    real(r8) :: grad(this%nelev)  ! computed gradient
    real(r8) :: grad_initial_guess(this%nelev)  ! initial guess for gradient
    real(r8) :: topo_interface(0:this%nelev)  ! elevations at interfaces between classes
    real(r8) :: h_lo(this%nelev) ! lower bounds for computing norms
    real(r8) :: h_hi(this%nelev) ! upper bounds for computing norms
    real(r8) :: dgrad(this%nelev) ! grad - grad_initial_guess
    real(r8) :: weight_grad(this%nelev) ! weight for dgrad in solution

    real(r8) :: diag(this%nelev-1) ! diagonal of tridiagonal matrix
    real(r8) :: subd(this%nelev-1) ! subdiagonal of tridiagonal matrix
    real(r8) :: supd(this%nelev-1) ! superdiagonal of tridiagonal matrix

    real(r8) :: b(this%nelev-1) ! rhs in A*x = b
    real(r8) :: A(this%nelev-1, this%nelev) ! matrix in A*x = b
    real(r8) :: A_AT(this%nelev-1, this%nelev-1) ! A * (transpose of A)
    real(r8) :: Tinv(this%nelev-1, this%nelev-1) ! inverse of tridiagonal matrix T = A * AT
    real(r8) :: AT_Tinv(this%nelev, this%nelev-1) ! (transpose of A) * Tinv
    real(r8) :: x_least_norm(this%nelev) ! least-norm solution x in A*x = b

    ! FIXME(wjs, 2016-04-27) Rename to nelev, or probably just delete
    integer :: n

    ! FIXME(wjs, 2016-04-26) Rename to ec
    integer :: i


    character(len=*), parameter :: subname = 'solve_for_vertical_gradients'
    !-----------------------------------------------------------------------

    field(:) = this%field(:,pt)
    topo(:) = this%topo(:,pt)

    ! FIXME(wjs, 2016-04-27) Rename to elevclass_bounds? Update: just delete this
    ! temporary variable.
    topo_interface(:) = this%elevclass_bounds(:)

    n = this%nelev

    ! FIXME(wjs, 2016-04-26) Extract method for the following two loops: returns
    ! weight_grad in each elevation class
    do i = 1, n

       if (i == 1) then
          ! If topo(1) is near the top of EC1, then the weight just includes twice the
          ! width of [topo(i) .. topo_interface(i)]
          h_lo(i) = max(topo_interface(i-1), (topo(i) - (topo_interface(i) - topo(i))))
       else
          h_lo(i) = topo_interface(i-1)
       end if

       if (i == n) then
          ! If topo(n) is near the bottom of EC N, then the weight just includes twice
          ! the width of [topo_interface(i-1) .. topo(i)]
          h_hi(i) = min(topo_interface(i), (topo(i) + (topo(i) - topo_interface(i-1))))
       else
          h_hi(i) = topo_interface(i)
       end if

    end do

    do i = 1, n
       ! set gradient weights based on h_hi - h_lo in each class
       weight_grad(i) = (h_hi(i) - h_lo(i)) / (h_hi(n) - h_lo(1))

       ! FIXME(wjs, 2016-04-26) Check that weight_grad > 0.
       !       If weight_grad is just slightly > 0, the matrix will be poorly conditioned.

    end do

    call this%calculator_initial_guess%get_gradients_one_point(pt, grad_initial_guess)

    !--------------------------------------------------------------------
    ! Set up matrix problem for gradient solution.
    ! The idea is to match field values at interfaces.
    !
    ! For each class n:
    !    field(n) + du(n)*grad(n) = field(n+1) - dl(n+1)*grad(n+1)
    !
    ! Rearrange to get
    !    du(n)*grad(n) + dl(n+1)*grad(n+1) = field(n+1) - field(n)
    ! 
    ! This is a bidiagonal matrix system:
    !
    !    | du(1)    dl(2)                        |   | grad(1) |    | field(2) - field(1) |
    !    |          du(2)   dl(3)                |   | grad(2) |    | field(3) - field(2) |
    !    |                  du(3)   dl(4)        | * | grad(3) | =  | field(4) - field(3) |
    !    |                          du(4)  dl(5) |   | grad(4) |    | field(5) - field(4) |
    !                                                | grad(5) |
    !
    ! The solution is underdetermined (4 equations, 5 unknowns).
    ! So we add an additional constraint: 
    ! Minimize the norm of the difference between the gradient in each class and the
    !  mean gradient, weighted by the range of the elevation class.
    ! That is, minimize the sum over i of (wt(i) * (grad(i) - grad_initial_guess(i)))^2.
    !
    ! The mean gradient is given by 
    !
    !                field(n) - field(1)
    !    grad_mean = __________________
    !                 topo(n) - topo(1)
    !
    ! The weights are
    !
    !         h_hi(i) - h_lo(i)
    ! wt(i) = _________________
    !         h_hi(n) - h_lo(1)
    !
    ! Putting in the weights and rearranging, we get
    !
    !         | wt(1) * dgrad(1) |    | field(2) - field(1) |       | wt(1) * grad_mean |
    !         | wt(2) * dgrad(2) |    | field(3) - field(2) |       | wt(2) * grad_mean |
    !     A * | wt(3) * dgrad(3) | =  | field(4) - field(3) | - A * | wt(3) * grad_mean |
    !         | wt(4) * dgrad(4) |    | field(5) - field(4) |       | wt(4) * grad_mean |
    !         | wt(5) * dgrad(5) |                                  | wt(5) * grad_mean |
    !
    ! where A is the bidiagonal matrix above, adjusted to include the gradient weights.
    ! E.g., du(1) becomes du(1)/wt(1), etc.
    ! This system is in the form A*x = b, where we want to solve for x.
    !
    ! It can be shown that the least-norm solution of A*x = b is given by
    !
    !    x = A^T (A*A^T)^{-1} b
    !
    ! So given A, we need to compute (A*A^T), takes its inverse, premultiply by A^T,
    ! and multiply the result by b.
    ! Given the least-norm solution, it is straightforward to compute the gradients.
    !--------------------------------------------------------------------

    ! Fill A
    ! A has (n-1) rows and n columns

    A(:,:) = 0._r8
    do i = 1, n-1
       A(i,i) = this%du(pt,i) / weight_grad(i)
       A(i,i+1) = this%dl(pt,i+1) / weight_grad(i+1)
    end do

    ! Compute A * A^T, a tridiagonal matrix of size (n-1)

    A_AT = matmul(A, transpose(A))

    ! Compute tridiagonal entries of (A * A^T)

    do i = 1, n-1
       diag(i) = A_AT(i,i)        ! diagonal
       if (i < n-1) then
          supd(i) = A_AT(i,i+1)   ! superdiagonal
          subd(i) = A_AT(i+1,i)   ! subdiagonal
       else
          supd(i) = 0._r8
          subd(i) = 0._r8
       end if
    end do

    ! Compute inverse of (A * A^T)
    ! Tinv has size (n-1,n-1)

    call tridiagonal_inverse(diag, supd, subd, Tinv)

    ! Premultiply by A^T
    ! A^T * Tinv has n rows and (n-1) columns

    AT_Tinv = matmul(transpose(A), Tinv)

    ! Compute the rhs vector b
    ! Start with the field differences

    do i = 1, n-1
       b(i) = field(i+1) - field(i)
    enddo

    ! Subtract (A * weight_grad*grad_initial_guess), a vector of size(n-1)

    b(:) = b(:) - matmul(A, weight_grad(:)*grad_initial_guess(:))

    ! Multiply AT_Tinv by b to get the least-norm solution
    ! b has size (n-1), x_least_norm has size n

    x_least_norm = matmul(AT_Tinv, b)

    ! Divide by the weighting factor to get dgrad
    dgrad(:) = x_least_norm(:) / weight_grad(:)

    ! Add dgrad to the target gradient to get the total gradient
    grad(:) = grad_initial_guess(:) + dgrad(:)

    ! Limit gradients
    call this%limit_gradients( &
         pt = pt, &
         grad_initial_guess = grad_initial_guess, &
         grad = grad)

    ! Finally, set class-level values
    this%vertical_gradient(:,pt) = grad(:)

  end subroutine solve_for_vertical_gradients

  !-----------------------------------------------------------------------
  subroutine limit_gradients(this, pt, grad_initial_guess, grad)
    !
    ! !DESCRIPTION:
    ! Limit the computed gradients for the given point
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_continuous_type), intent(inout) :: this
    integer, intent(in) :: pt
    
    ! All of these arguments should have size this%nelev

    ! we'll back off to these initial guesses if there is a problem with grad in a given
    ! elevation class
    real(r8), intent(in) :: grad_initial_guess(:)

    ! upon input, grad contains the current gradient estimates; upon output, it is
    ! modified to be limited
    real(r8), intent(inout) :: grad(:)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: field(this%nelev) ! mean field value of each elevation class
    integer :: ec
    real(r8) :: val_at_lb  ! value at the lower bound interface of an elevation class
    real(r8) :: val_at_ub  ! value at the upper bound interface of an elevation class
    logical :: val_at_lb_outside_bounds
    logical :: val_at_ub_outside_bounds

    character(len=*), parameter :: subname = 'limit_gradients'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(grad_initial_guess) == (/this%nelev/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(grad) == (/this%nelev/)), errMsg(__FILE__, __LINE__))

    field(:) = this%field(:,pt)

    ! Set gradient to 0 in lowest and highest elevation class
    grad(1) = 0._r8
    this%limited_to_zero(1, pt) = .true.
    grad(this%nelev) = 0._r8
    this%limited_to_zero(this%nelev, pt) = .true.

    do ec = 2, this%nelev - 1
       val_at_lb = field(ec) - (this%dl(pt,ec) * grad(ec))
       val_at_lb_outside_bounds = is_outside_bounds(val_at_lb, field(ec), field(ec-1))

       val_at_ub = field(ec) + (this%du(pt,ec) * grad(ec))
       val_at_ub_outside_bounds = is_outside_bounds(val_at_ub, field(ec), field(ec+1))

       if (val_at_lb_outside_bounds .or. val_at_ub_outside_bounds) then
          grad(ec) = grad_initial_guess(ec)
          this%limited_to_initial_guess(ec, pt) = .true.
       end if
    end do

  contains
    pure logical function is_outside_bounds(val, bound1, bound2)
      ! Returns true if val is outside the interval given by bound1 and bound2
      real(r8), intent(in) :: val, bound1, bound2

      if (val < bound1 .and. val < bound2) then
         is_outside_bounds = .true.
      else if (val > bound1 .and. val > bound2) then
         is_outside_bounds = .true.
      else
         is_outside_bounds = .false.
      end if
    end function is_outside_bounds
  end subroutine limit_gradients

  !-----------------------------------------------------------------------
  function dl(this, pt, ec)
    !
    ! !DESCRIPTION:
    ! Return lower half-width of elevation class ec in point pt
    !
    ! !ARGUMENTS:
    real(r8) :: dl  ! function result
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    integer, intent(in) :: pt
    integer, intent(in) :: ec
    !-----------------------------------------------------------------------

    dl = this%topo(ec, pt) - this%elevclass_bounds(ec - 1)

  end function dl


  !-----------------------------------------------------------------------
  function du(this, pt, ec)
    !
    ! !DESCRIPTION:
    ! Return upper half-width of elevation class ec in point pt
    !
    ! !ARGUMENTS:
    real(r8) :: du  ! function result
    class(vertical_gradient_calculator_continuous_type), intent(in) :: this
    integer, intent(in) :: pt
    integer, intent(in) :: ec
    !-----------------------------------------------------------------------

    du = this%elevclass_bounds(ec) - this%topo(ec, pt)

  end function du


end module vertical_gradient_calculator_continuous
