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
  
  implicit none
  private

  public :: vertical_gradient_calculator_constant_type

  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_constant_type
     private
     integer :: num_points
     real(r8) :: gradient
   contains
     procedure :: calc_vertical_gradient
  end type vertical_gradient_calculator_constant_type

  interface vertical_gradient_calculator_constant_type
     module procedure constructor
  end interface vertical_gradient_calculator_constant_type

contains

  !-----------------------------------------------------------------------
  function constructor(num_points, gradient) result(this)
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
    real(r8), intent(in) :: gradient
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%num_points = num_points
    this%gradient = gradient
    
  end function constructor

  !-----------------------------------------------------------------------
  subroutine calc_vertical_gradient(this, elevation_class, vertical_gradient)
    !
    ! !DESCRIPTION:
    ! Calculate the vertical gradient for all points
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_constant_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! vertical_gradient should already be allocated to the appropriate size
    real(r8), intent(out) :: vertical_gradient(:)
    !
    ! !LOCAL VARIABLES:
    integer :: grid_cell
    
    character(len=*), parameter :: subname = 'calc_vertical_gradient'
    !-----------------------------------------------------------------------

    SHR_ASSERT((size(vertical_gradient) == this%num_points), subname//': wrong size for vertical gradient')

    do grid_cell = 1, this%num_points
       vertical_gradient(grid_cell) = this%gradient * elevation_class * grid_cell**2
    end do
  end subroutine calc_vertical_gradient

end module vertical_gradient_calculator_constant
