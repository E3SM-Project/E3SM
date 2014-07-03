module shr_sys_mod

! This is a fake module for testing debug_utils.
! It contains only a flush, and a fake abort
! method.

use shr_kind_mod, only: &
     shr_kind_in

implicit none
private
save

! Fake abort
public :: shr_sys_abort

! Test if abort was called, then reset flag.
public :: pull_aborted

! Still want a real flush available.
public :: shr_sys_flush

! Stores whether shr_sys_abort was called.
logical :: aborted = .false.

contains

subroutine shr_sys_abort(string, rc)

  use shr_kind_mod, only: shr_kind_in

  character(*), optional :: string
  integer(shr_kind_in), optional :: rc

  aborted = .true.

end subroutine shr_sys_abort

function pull_aborted() result(flag_out)
  logical :: flag_out

  flag_out = aborted
  aborted = .false.

end function pull_aborted

SUBROUTINE shr_sys_flush(unit)

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

#ifndef CPRLAHEY
   flush(unit)
#else
   call flush(unit)
#endif

END SUBROUTINE shr_sys_flush

end module shr_sys_mod
