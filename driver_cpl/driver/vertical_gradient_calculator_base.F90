module vertical_gradient_calculator_base

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines an abstract base class for computing the vertical gradient of a
  ! field.
  !
  ! Usage:
  !
  ! - First call calc_gradients
  !
  ! - Then can query the computed vertical gradients using the other methods

  use seq_comm_mct, only : logunit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  
  implicit none
  private

  public :: vertical_gradient_calculator_base_type

  type, abstract :: vertical_gradient_calculator_base_type
   contains
     ! Calculate the vertical gradients for all points and all elevation classes
     procedure(calc_gradients_interface), deferred :: calc_gradients

     ! Get the vertical gradients for all points for a single elevation class
     procedure(get_gradients_one_class_interface), deferred :: get_gradients_one_class

     ! Get the vertical gradients for all elevation classes for a single point
     procedure(get_gradients_one_point_interface), deferred :: get_gradients_one_point

     ! These routines are utility methods for derived classes; they should not be called
     ! by clients of this class.
     procedure, nopass :: check_elevclass_bounds_monotonic_increasing


  end type vertical_gradient_calculator_base_type

  abstract interface
     subroutine calc_gradients_interface(this)
       import :: vertical_gradient_calculator_base_type
       class(vertical_gradient_calculator_base_type), intent(inout) :: this
     end subroutine calc_gradients_interface

     subroutine get_gradients_one_class_interface(this, elevation_class, gradients)
       import :: vertical_gradient_calculator_base_type
       import :: r8
       class(vertical_gradient_calculator_base_type), intent(in) :: this
       integer, intent(in) :: elevation_class

       ! vertical_gradient should already be allocated to the appropriate size
       real(r8), intent(out) :: gradients(:)
     end subroutine get_gradients_one_class_interface

     subroutine get_gradients_one_point_interface(this, point, gradients)
       import :: vertical_gradient_calculator_base_type
       import :: r8
       class(vertical_gradient_calculator_base_type), intent(in) :: this
       integer, intent(in) :: point

       ! vertical_gradient should already be allocated to the appropriate size
       real(r8), intent(out) :: gradients(:)
     end subroutine get_gradients_one_point_interface
  end interface

contains

  !-----------------------------------------------------------------------
  subroutine check_elevclass_bounds_monotonic_increasing(elevclass_bounds)
    !
    ! !DESCRIPTION:
    ! Ensure that elevclass_bounds are monotonically increasing; abort if there is a
    ! problem
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: elevclass_bounds(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'check_elevclass_bounds'
    !-----------------------------------------------------------------------

    do i = 2, size(elevclass_bounds)
       if (elevclass_bounds(i-1) >= elevclass_bounds(i)) then
          write(logunit,*) subname, ': ERROR: elevclass_bounds must be monotonically increasing'
          write(logunit,*) 'elevclass_bounds = ', elevclass_bounds
          call shr_sys_abort(subname//': ERROR: elevclass_bounds must be monotonically increasing')
       end if
    end do

  end subroutine check_elevclass_bounds_monotonic_increasing



end module vertical_gradient_calculator_base
