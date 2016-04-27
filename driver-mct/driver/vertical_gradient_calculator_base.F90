module vertical_gradient_calculator_base

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines an abstract base class for computing the vertical gradient of a
  ! field.

  use seq_comm_mct, only : logunit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod, only : shr_sys_abort
  
  implicit none
  private
  
  public :: vertical_gradient_calculator_base_type
  
  type, abstract :: vertical_gradient_calculator_base_type
   contains
     ! Calculate the vertical gradient for all points, for a given elevation class
     procedure(calc_vertical_gradient_interface), deferred :: calc_vertical_gradient

     ! These routines are utility methods for derived classes; they should not be called
     ! by clients of this class.
     procedure, nopass :: check_elevclass_bounds_monotonic_increasing


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
