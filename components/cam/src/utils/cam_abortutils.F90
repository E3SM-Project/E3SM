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
  use shr_kind_mod, only: shr_kind_in
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

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!========================================================================

   subroutine endrun(string,rc)
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
    character(len=*)    , intent(in), optional :: string ! error message 
                                                         ! string
    integer(shr_kind_in), intent(in), optional :: rc     ! error code
    !
    ! Local workspace
    !
    character(len=*), parameter :: &
      cam_string = "See atm.log for error description" ! generic string
    !--------------------------------------------------------------------

    if (present(string)) then
      if (present(rc)) then
        call shr_sys_abort(string,rc)
      else
        call shr_sys_abort(string)
      endif
    else
      if (present(rc)) then
        call shr_sys_abort(cam_string,rc)
      else
        call shr_sys_abort(cam_string)
      endif
    endif

   end subroutine

end module cam_abortutils
