module vertical_gradient_calculator_continuousNoLimit
  ! This module provides a type that inherits from
  ! vertical_gradient_calculator_continuous, overriding the limiting to result in no
  ! limiting of the initially-computed gradients.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use vertical_gradient_calculator_base, only : &
       vertical_gradient_calculator_base_type
  use vertical_gradient_calculator_continuous, only : &
       vertical_gradient_calculator_continuous_type

  implicit none
  private

  public :: vgc_continuousNoLimit_type

  type, extends(vertical_gradient_calculator_continuous_type) :: &
       vgc_continuousNoLimit_type
     private

   contains
     procedure :: limit_gradients
  end type vgc_continuousNoLimit_type

  interface vgc_continuousNoLimit_type
     module procedure constructor
  end interface vgc_continuousNoLimit_type

contains

  function constructor(field, topo, elevclass_bounds, calculator_initial_guess) result(this)
    type(vgc_continuousNoLimit_type) :: this
    real(r8), intent(in) :: field(:,:)  ! field(i,j) is point i, elevation class j
    real(r8), intent(in) :: topo(:,:)   ! topo(i,j) is point i, elevation class j
    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8)         , intent(in) :: elevclass_bounds(0:)
    class(vertical_gradient_calculator_base_type), intent(in) :: calculator_initial_guess

    this%vertical_gradient_calculator_continuous_type = &
         vertical_gradient_calculator_continuous_type( &
         field, topo, elevclass_bounds, calculator_initial_guess)
  end function constructor

  subroutine limit_gradients(this, pt, grad_initial_guess, grad)
    class(vgc_continuousNoLimit_type), intent(in) :: this
    integer, intent(in) :: pt
    real(r8), intent(in) :: grad_initial_guess(:)
    real(r8), intent(inout) :: grad(:)

    ! do nothing
  end subroutine limit_gradients

end module vertical_gradient_calculator_continuousNoLimit
