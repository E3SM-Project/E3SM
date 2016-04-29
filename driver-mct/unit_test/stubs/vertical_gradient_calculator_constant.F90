module vertical_gradient_calculator_constant

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type that is
  ! useful for unit testing. It computes the gradient as a constant times the elevation
  ! class index times (the grid cell index squared).

#include "shr_assert.h"
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod, only : shr_sys_abort
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  
  implicit none
  private

  public :: vertical_gradient_calculator_constant_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_constant_type
     private
     integer :: num_points
     integer :: nelev
     real(r8) :: gradient
     real(r8), allocatable :: vertical_gradient(:,:)  ! [point, elev classs]
     logical :: calculated
   contains
     procedure :: calc_gradients
     procedure :: get_gradients_one_class
  end type vertical_gradient_calculator_constant_type

  interface vertical_gradient_calculator_constant_type
     module procedure constructor
  end interface vertical_gradient_calculator_constant_type

contains

  !-----------------------------------------------------------------------
  function constructor(num_points, nelev, gradient) result(this)
    !
    ! !DESCRIPTION:
    ! Create a new vertical_gradient_calculator_constant_type object.
    !
    ! The returned gradient will be (gradient) * (elevation_class) * (grid_cell_index)^2
    !
    ! num_points gives the number of points for which a gradient is needed (e.g., if
    ! computing the vertical gradient on the land domain, then num_points is the number
    ! of land points).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_constant_type) :: this  ! function result
    integer, intent(in) :: num_points
    integer, intent(in) :: nelev
    real(r8), intent(in) :: gradient
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%calculated = .false.
    this%num_points = num_points
    this%nelev = nelev
    this%gradient = gradient

    allocate(this%vertical_gradient(num_points, nelev))
    this%vertical_gradient(:,:) = nan

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
    class(vertical_gradient_calculator_constant_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: pt, ec

    character(len=*), parameter :: subname = 'calc_gradients'
    !-----------------------------------------------------------------------

    SHR_ASSERT(.not. this%calculated, 'gradients already calculated')

    do ec = 1, this%nelev
       do pt = 1, this%num_points
          this%vertical_gradient(pt, ec) = this%gradient * ec * pt**2
       end do
    end do

    this%calculated = .true.

  end subroutine calc_gradients


  !-----------------------------------------------------------------------
  subroutine get_gradients_one_class(this, elevation_class, gradients)
    !
    ! !DESCRIPTION:
    ! Calculate the vertical gradient for all points
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_constant_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! gradients should already be allocated to the appropriate size
    real(r8), intent(out) :: gradients(:)
    !
    ! !LOCAL VARIABLES:
    integer :: grid_cell
    
    character(len=*), parameter :: subname = 'get_gradients_one_class'
    !-----------------------------------------------------------------------

    SHR_ASSERT(this%calculated, 'gradients not yet calculated')
    SHR_ASSERT(elevation_class <= this%nelev, subname//': elevation class exceeds bounds')
    SHR_ASSERT((size(gradients) == this%num_points), subname//': wrong size for vertical gradient')

    gradients(:) = this%vertical_gradient(:, elevation_class)
  end subroutine get_gradients_one_class

end module vertical_gradient_calculator_constant
