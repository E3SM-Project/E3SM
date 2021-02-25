module shr_abort_mod

  ! This is a replacement for shr_abort_mod that throws a pfunit exception rather than
  ! aborting

  use shr_kind_mod, only : shr_kind_in
  use pfunit_mod, only : throw

  implicit none
  private

  public :: shr_abort_abort ! Replacement for shr_abort_abort that throws a pfunit exception rather than aborting

  public :: shr_abort_backtrace ! Just to satisfy the public interface of shr_abort_abort

contains

  subroutine shr_abort_abort(string,rc)
    ! Replacement for shr_abort_abort that throws a pfunit exception rather than aborting
    !
    ! This can be used to test expected errors (i.e., failure testing).
    !
    ! If this occurs within a pFUnit-based test:
    !
    ! - If you have code like:
    !
    !   @assertExceptionRaised(expected_message)
    !
    !   then your test will pass if the actual message in the 'throw' call (including the
    !   'ABORTED: ' prefix) matches expected_message; it will fail if the actual message
    !   doesn't match the expected message
    !
    ! - If you don't have
    !
    !   @assertExceptionRaised
    !
    !   or
    !
    !   call assertExceptionRaised
    !
    !   then this will result in the given pFUnit test failing.

    !----- arguments -----
    character(len=*)    , intent(in), optional :: string  ! error message string
    integer(shr_kind_in), intent(in), optional :: rc      ! error code

    !----- locals -----
    integer(shr_kind_in) :: my_rc

    ! Prevent compiler spam about unused variables.
    if (.false.) my_rc = rc

    call throw("ABORTED: "//trim(string))
  end subroutine shr_abort_abort

  subroutine shr_abort_backtrace()
    ! Just to satisfy the public interface of shr_abort_abort
    !
    ! Does nothing

  end subroutine shr_abort_backtrace

end module shr_abort_mod
