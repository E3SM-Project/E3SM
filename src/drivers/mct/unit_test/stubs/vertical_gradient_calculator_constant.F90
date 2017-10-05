module vertical_gradient_calculator_specified

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type that is
  ! useful for unit testing. It returns a specified vertical gradient.
  !
  ! This module also provides convenience functions for creating a
  ! vertical_gradient_calculator_specified_type object with various functional forms.

  ! It computes the gradient as a constant times the elevation
  ! class index times (the grid cell index squared).

#include "shr_assert.h"
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod, only : shr_sys_abort
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  
  implicit none
  private

  public :: vertical_gradient_calculator_specified_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_specified_type
     private
     integer :: num_points
     integer :: nelev
     real(r8), allocatable :: vertical_gradient(:,:)  ! [point, elev classs]
     logical :: calculated
   contains
     procedure :: calc_gradients
     procedure :: get_gradients_one_class
     procedure :: get_gradients_one_point
  end type vertical_gradient_calculator_specified_type

  interface vertical_gradient_calculator_specified_type
     module procedure constructor
  end interface vertical_gradient_calculator_specified_type

  ! Creates a calculator where the gradient in ec i, pt j is gradient * i * j^2
  public :: vgc_specified_ec_times_ptSquared

  ! Creates a calculator where the gradient is constant for each point, set as the mean
  ! slope from the lowest to highest elev class
  public :: vgc_specified_mean_slope
contains

  !-----------------------------------------------------------------------
  function vgc_specified_ec_times_ptSquared(num_points, nelev, gradient) &
       result(calculator)
    !
    ! !DESCRIPTION:
    ! Creates a calculator where the gradient in ec i, pt j is gradient * i * j^2
    !
    ! num_points gives the number of points for which a gradient is needed (e.g., if
    ! computing the vertical gradient on the land domain, then num_points is the number
    ! of land points).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_specified_type) :: calculator  ! function result
    integer, intent(in) :: num_points
    integer, intent(in) :: nelev
    real(r8), intent(in) :: gradient
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: gradients(:,:)
    integer :: pt, ec

    character(len=*), parameter :: subname = 'vgc_specified_ec_times_ptSquared'
    !-----------------------------------------------------------------------

    allocate(gradients(num_points, nelev))

    do ec = 1, nelev
       do pt = 1, num_points
          gradients(pt, ec) = gradient * ec * pt**2
       end do
    end do

    calculator = vertical_gradient_calculator_specified_type(gradients)

  end function vgc_specified_ec_times_ptSquared

  !-----------------------------------------------------------------------
  function vgc_specified_mean_slope(data, topo) result(calculator)
    !
    ! !DESCRIPTION:
    ! Creates a calculator where the gradient is constant for all elevation classes -
    ! though can differ for each point. Specifically, it is set to the mean slope from
    ! the lowest to highest elev class
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_specified_type) :: calculator  ! function result
    real(r8), intent(in) :: data(:,:)  ! [pt, ec]
    real(r8), intent(in) :: topo(:,:)  ! [pt, ec]
    !
    ! !LOCAL VARIABLES:
    integer :: num_points
    integer :: nelev
    real(r8), allocatable :: gradients(:,:)
    integer pt

    character(len=*), parameter :: subname = 'vgc_specified_mean_slope'
    !-----------------------------------------------------------------------

    num_points = size(data,1)
    nelev = size(data,2)
    SHR_ASSERT_ALL((ubound(topo) == (/num_points, nelev/)), 'bad size for topo')

    allocate(gradients(num_points, nelev))

    do pt = 1, num_points
       gradients(pt, :) = (data(pt,nelev) - data(pt,1)) / &
            (topo(pt,nelev) - topo(pt,1))
    end do

    calculator = vertical_gradient_calculator_specified_type(gradients)

  end function vgc_specified_mean_slope

  !-----------------------------------------------------------------------
  function constructor(gradients) result(this)
    !
    ! !DESCRIPTION:
    ! Create a new vertical_gradient_calculator_specified_type object.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_specified_type) :: this  ! function result
    real(r8), intent(in) :: gradients(:,:)  ! [pt, ec]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%calculated = .false.
    this%num_points = size(gradients, 1)
    this%nelev = size(gradients, 2)

    allocate(this%vertical_gradient(this%num_points, this%nelev))
    this%vertical_gradient(:,:) = gradients(:,:)

  end function constructor

  !-----------------------------------------------------------------------
  subroutine calc_gradients(this)
    !
    ! !DESCRIPTION:
    ! Calculate the vertical gradients
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_specified_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'calc_gradients'
    !-----------------------------------------------------------------------

    SHR_ASSERT(.not. this%calculated, 'gradients already calculated')

    ! Nothing to do in this stub

    this%calculated = .true.

  end subroutine calc_gradients


  !-----------------------------------------------------------------------
  subroutine get_gradients_one_class(this, elevation_class, gradients)
    !
    ! !DESCRIPTION:
    ! Return the vertical gradients for one elevation class, for all points
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_specified_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'get_gradients_one_class'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, 'gradients not yet calculated')
    SHR_ASSERT(elevation_class <= this%nelev, subname//': elevation class exceeds bounds')
    SHR_ASSERT((size(gradients) == this%num_points), subname//': wrong size for vertical gradient')

    gradients(:) = this%vertical_gradient(:, elevation_class)
  end subroutine get_gradients_one_class

  !-----------------------------------------------------------------------
  subroutine get_gradients_one_point(this, point, gradients)
    !
    ! !DESCRIPTION:
    ! Return the vertical gradient for all elevation classes, for one point
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_specified_type), intent(in) :: this
    integer, intent(in) :: point

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_gradients_one_class'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, 'gradients not yet calculated')
    SHR_ASSERT(point <= this%num_points, subname//': elevation class exceeds bounds')
    SHR_ASSERT((size(gradients) == this%nelev), subname//': wrong size for vertical gradient')

    gradients(:) = this%vertical_gradient(point, :)
  end subroutine get_gradients_one_point

end module vertical_gradient_calculator_specified
