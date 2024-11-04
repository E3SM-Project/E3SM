module cam_abortutils

!----------------------------------------------------------------------- 
! 
! Purpose: This module provides CAM an interface to the model
!          stopping mechanisms.
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
  use shr_kind_mod, only: shr_kind_in, SHR_KIND_CL
  use shr_sys_mod,  only: shr_sys_abort

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
  implicit none
  private
  save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
  public endrun
  public handle_allocate_error

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!========================================================================

   subroutine endrun(string)
    !--------------------------------------------------------------------
    ! 
    ! Purpose: CAM interface to consistent stopping mechanism
    ! 
    ! Method: Calls shr_sys_abort
    ! 
    !--------------------------------------------------------------------

    implicit none

    !
    ! Arguments
    !
    character(len=*), intent(in), optional :: string   ! error message

    !
    ! Local workspace
    !
    character(len=*), parameter :: &
      cam_string = "See atm.log for error description" ! generic message
    !--------------------------------------------------------------------

    if (present(string)) then
      call shr_sys_abort(string)
    else
      call shr_sys_abort(cam_string)
    endif

   end subroutine

   subroutine handle_allocate_error(retval, subname, fieldname)
      ! if <retval> is not zero, generate an error message and abort
      ! Dummy arguments
      integer,          intent(in) :: retval
      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: fieldname
      ! Local variable
      character(len=SHR_KIND_CL)   :: errmsg

      if (retval /= 0) then
         write(errmsg, '(4a,i0)') trim(subname), ' error allocating ',        &
              trim(fieldname), ', error = ', retval
         call endrun(errmsg)
      end if
   end subroutine handle_allocate_error

end module cam_abortutils
