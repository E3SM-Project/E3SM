module vertical_gradient_calculator_2nd_order

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type for
  ! computing vertical gradients using a second-order centered difference.
  !
  ! If the topo values are nearly equal across the gradient (i.e., denominator is near 0),
  ! returns a gradient of 0.
  !
  ! If there is only one elevation class, returns a gradient of 0.

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_sys_mod, only : shr_sys_abort
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  private

  public :: vertical_gradient_calculator_2nd_order_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_2nd_order_type
     private

     integer :: nelev  ! number of elevation classes
     integer :: num_points
     real(r8), allocatable :: field(:,:)   ! field(i,j) is point i, elevation class j
     real(r8), allocatable :: topo(:,:)    ! topo(i,j) is point i, elevation class j

     ! Bounds of each elevation class. This array has one more element than the number of
     ! elevation classes, since it contains lower and upper bounds for each elevation
     ! class. The indices go 0:nelev. These bounds
     ! are guaranteed to be monotonically increasing.
     real(r8), allocatable :: elevclass_bounds(:)

     ! precomputed vertical gradients; vertical_gradient(i,j) is point i, elevation class
     ! j
     real(r8), allocatable :: vertical_gradient(:,:)

     logical :: calculated  ! whether gradients have been calculated yet

   contains
     procedure :: calc_gradients
     procedure :: get_gradients_one_class
     procedure :: get_gradients_one_point

     procedure, private :: check_topo ! check topographic heights
     procedure, private :: limit_gradient

  end type vertical_gradient_calculator_2nd_order_type

  interface vertical_gradient_calculator_2nd_order_type
     module procedure constructor
  end interface vertical_gradient_calculator_2nd_order_type

contains

  !-----------------------------------------------------------------------
  function constructor(field, topo, elevclass_bounds) result(this)
    !
    ! !DESCRIPTION:
    ! Creates a vertical_gradient_calculator_2nd_order_type object.
    !
    ! Pre-condition: elevclass_bounds must be monotonically increasing.
    !
    ! Pre-condition: Topographic heights must all lie inside the bounds of their
    ! respective elevation class (given by elevclass_bounds), with the possible exception
    ! of the lowest elevation class (topographic heights can lie below the arbitrary lower
    ! bound of the elevation class) and the highest elevation class (topographic heights
    ! can lie above the arbitrary upper bound of the elevation class). (This pre-condition
    ! is mainly important for the sake of calculating the limiter.)
    !
    ! TODO(wjs, 2016-04-21) Currently the topographic heights pre-condition is not
    ! checked: see below.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_2nd_order_type) :: this  ! function result
    real(r8), intent(in) :: field(:,:)  ! field(i,j) is point i, elevation class j
    real(r8), intent(in) :: topo(:,:)   ! topo(i,j) is point i, elevation class j

    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8), intent(in) :: elevclass_bounds(0:)
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

    ! (In principle, we could also handle monotonically decreasing elevclass_bounds, but
    ! that would require generalizing some code, such as in check_topo.)
    call this%check_elevclass_bounds_monotonic_increasing(this%elevclass_bounds)

    allocate(this%field(this%num_points, this%nelev))
    this%field(:,:) = field(:,:)
    allocate(this%topo(this%num_points, this%nelev))
    this%topo(:,:) = topo(:,:)

    ! TODO(wjs, 2016-04-21) Uncomment this call to check_topo. It's important for
    ! topographic heights to be within bounds in order for the limiter to be applied
    ! correctly. However, this currently isn't the case for some of the old TG forcing
    ! data. At a glance, it looks like the problems are just outside of Greenland, so this
    ! should be okay. When we have new TG forcing data, we should try uncommenting this
    ! call to check_topo.
    !
    ! Alternatively, we could change check_topo to set a flag for each point saying
    ! whether topo values are bad for that point. Then, when computing gradients, set
    ! them to 0 for all points with bad topo values. (We did something similar for the
    ! now-deleted vertical_gradient_calculator_continuous.) However, longer-term, an
    ! abort may be more appropriate rather than silently setting gradients to 0.

    ! call this%check_topo()

    allocate(this%vertical_gradient(this%num_points, this%nelev))
    this%vertical_gradient(:,:) = nan

  end function constructor

  !-----------------------------------------------------------------------
  subroutine calc_gradients(this)
    !
    ! !DESCRIPTION:
    ! Calculates all vertical gradients
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    ! Tolerance for considering two topo values to be nearly equal
    real(r8), parameter :: topo_equality_tolerance = 1.e-13_r8

    integer :: i
    integer :: elevation_class
    integer :: ec_low   ! elevation class index to use as the lower bound of the gradient
    integer :: ec_high  ! elevation class index to use as the upper bound of the gradient
    logical :: two_sided ! true if we're estimating the gradient with a two-sided difference

    character(len=*), parameter :: subname = 'calc_gradients'
    !-----------------------------------------------------------------------

    if (this%calculated) then
       ! nothing to do
       return
    end if

    if (this%nelev == 1) then
       this%vertical_gradient(:,:) = 0._r8
       two_sided = .false.

    else

       do elevation_class = 1, this%nelev
          if (elevation_class == 1) then
             ec_low = elevation_class
             ec_high = elevation_class + 1
             two_sided = .false.
          else if (elevation_class == this%nelev) then
             ec_low = elevation_class - 1
             ec_high = elevation_class
             two_sided = .false.
          else
             ec_low = elevation_class - 1
             ec_high = elevation_class + 1
             two_sided = .true.
          end if

          do i = 1, this%num_points
             if (abs(this%topo(i, ec_high) - this%topo(i, ec_low)) < topo_equality_tolerance) then
                this%vertical_gradient(i, elevation_class) = 0._r8
             else
                this%vertical_gradient(i, elevation_class) = &
                     (this%field(i, ec_high) - this%field(i, ec_low)) / &
                     (this%topo (i, ec_high) - this%topo (i, ec_low))
             end if
          end do

          if (two_sided) then
             call this%limit_gradient(elevation_class, ec_low, ec_high, &
                  this%vertical_gradient(:,elevation_class))
          end if
       end do

    end if
    
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
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_gradients_one_class'
    !-----------------------------------------------------------------------

    ! Assert pre-conditions

    SHR_ASSERT(this%calculated, errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(gradients) == this%num_points), errMsg(__FILE__, __LINE__))

    if (elevation_class < 1 .or. elevation_class > this%nelev) then
       write(logunit,*) subname, ': ERROR: elevation class out of bounds: ', &
            elevation_class, this%nelev
       call shr_sys_abort(subname//': ERROR: elevation class out of bounds')
    end if

    gradients(:) = this%vertical_gradient(:, elevation_class)

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
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
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

    gradients(:) = this%vertical_gradient(point, :)

  end subroutine get_gradients_one_point

  !-----------------------------------------------------------------------
  subroutine check_topo(this)
    !
    ! !DESCRIPTION:
    ! Check topographic heights; abort if there is a problem
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
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
    integer :: i

    ! Absolute tolerance for error checks. This is chosen so that it allows for
    ! double-precision roundoff-level errors on values of order 10,000.
    real(r8), parameter :: tol = 1.e-10_r8

    character(len=*), parameter :: subname = 'check_topo'
    !-----------------------------------------------------------------------

    do elevclass = 1, this%nelev
       if (elevclass > 1) then
          do i = 1, this%num_points
             if (this%topo(i,elevclass) - this%elevclass_bounds(elevclass-1) < -tol) then
                write(logunit,*) subname, ': ERROR: topo lower than lower bound of elevation class:'
                write(logunit,*) 'i, elevclass, topo, lower_bound = ', &
                     i, elevclass, this%topo(i,elevclass), this%elevclass_bounds(elevclass-1)
                call shr_sys_abort(subname//': ERROR: topo lower than lower bound of elevation class')
             end if
          end do
       end if

       if (elevclass < this%nelev) then
          do i = 1, this%num_points
             if (this%topo(i,elevclass) - this%elevclass_bounds(elevclass) > tol) then
                write(logunit,*) subname, ': ERROR: topo higher than upper bound of elevation class:'
                write(logunit,*) 'i, elevclass, topo, upper_bound = ', &
                     i, elevclass, this%topo(i,elevclass), this%elevclass_bounds(elevclass)
                call shr_sys_abort(subname//': ERROR: topo higher than upper bound of elevation class')
             end if
          end do
       end if
    end do

  end subroutine check_topo

  !-----------------------------------------------------------------------
  subroutine limit_gradient(this, k, ec_low, ec_high, vertical_gradient)
    !
    ! !DESCRIPTION:
    ! Limit the gradient: Ensure that the interface values lie inside the range defined
    ! by the max and min of the mean values in this class and its 2 adjacent neighbors.
    !
    ! Should only be called for two-sided differences (ec_low < k < ec_high)
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    integer , intent(in) :: k       ! elevation class index
    integer , intent(in) :: ec_low  ! elevation class index used as the lower bound of the gradient
    integer , intent(in) :: ec_high ! elevation class index used as the upper bound of the gradient
    real(r8), intent(inout) :: vertical_gradient(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    real(r8) :: deviation_high
    real(r8) :: deviation_low
    real(r8) :: deviation_max
    real(r8) :: deviation_min
    real(r8) :: diff_max
    real(r8) :: diff_min
    real(r8) :: factor1
    real(r8) :: factor2
    real(r8) :: limiting_factor

    character(len=*), parameter :: subname = 'limit_gradient'
    !-----------------------------------------------------------------------

    ! Basic idea: In 1D with a linear reconstruction, the extreme values of the data will
    ! lie at the interfaces between adjacent elevation classes. The interface values
    ! should not lie outside the range defined by the max and min of the mean values in
    ! this class and its 2 adjacent neighbors.

    ! This code only works correctly if we're doing a two-sided difference (otherwise,
    ! one of diff_min or diff_max will be 0, leading to 0 gradient - when in fact we
    ! don't want to do any limiting for a one-sided difference).
    SHR_ASSERT(ec_low < k, subname//': Only works for two-sided difference: must have ec_low < k')
    SHR_ASSERT(ec_high > k, subname//': Only works for two-sided difference: must have ec_high > k')

    do i = 1, this%num_points
       ! First compute the max and min values of the deviation of the data from its mean
       ! value. With a linear gradient, the max differences must lie at the adjacent
       ! interfaces.
       deviation_high = vertical_gradient(i) * (this%elevclass_bounds(k) - this%topo(i,k))
       deviation_low = vertical_gradient(i) *  (this%elevclass_bounds(k-1) - this%topo(i,k))
       deviation_max = max(deviation_high, deviation_low)
       deviation_min = min(deviation_high, deviation_low)

       ! Now compute the max and min of the data in the cell and its nearest neighbors.
       ! (Actually, the difference between this max/min value and the mean value in the
       ! current class.)
       diff_max = max(this%field(i,ec_high), this%field(i,k), this%field(i,ec_low)) - this%field(i,k)
       diff_min = min(this%field(i,ec_high), this%field(i,k), this%field(i,ec_low)) - this%field(i,k)

       ! Now limit the gradient using the information computed above.

       if (abs(deviation_min) > 0._r8) then
          factor1 = max(0._r8, diff_min/deviation_min)
       else
          factor1 = 1._r8
       endif
       
       if (abs(deviation_max) > 0._r8) then
          factor2 = max(0._r8, diff_max/deviation_max)
       else
          factor2 = 1._r8
       endif

       ! limiting factor will lie between 0 and 1
       limiting_factor = min(1._r8, factor1, factor2)
       vertical_gradient(i) = vertical_gradient(i) * limiting_factor
    end do

  end subroutine limit_gradient

end module vertical_gradient_calculator_2nd_order

