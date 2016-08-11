module vertical_gradient_calculator_base

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines an abstract base class for computing the vertical gradient of a
  ! field.

  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none
  private
  
  public :: vertical_gradient_calculator_base_type
  
  type, abstract :: vertical_gradient_calculator_base_type
   contains
     ! Calculate the vertical gradient for all points, for a given elevation class
     procedure(calc_vertical_gradient_interface), deferred :: calc_vertical_gradient
  end type vertical_gradient_calculator_base_type
  
  abstract interface
     subroutine calc_vertical_gradient_interface(this, elevation_class, vertical_gradient)
       import :: vertical_gradient_calculator_base_type
       import :: r8
       class(vertical_gradient_calculator_base_type), intent(in) :: this
       integer, intent(in) :: elevation_class

       ! vertical_gradient should already be allocated to the appropriate size
       real(r8), intent(out) :: vertical_gradient(:)
     end subroutine calc_vertical_gradient_interface
  end interface

end module vertical_gradient_calculator_base
