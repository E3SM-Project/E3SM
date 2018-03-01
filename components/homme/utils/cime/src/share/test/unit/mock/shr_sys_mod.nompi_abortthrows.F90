module shr_sys_mod

! This is a mock version of shr_sys_mod.
! It contains only a few routines that are needed, and an abort method that throws a pFUnit
! exception instead of actually aborting.

use shr_kind_mod, only: &
     shr_kind_in, shr_kind_r8

use shr_abort_mod, only : shr_sys_abort => shr_abort_abort

implicit none
private
save

! Fake abort
! Imported from shr_abort_mod and republished with renames. Other code that wishes to use
! these routines should use these shr_sys names rather than directly using the routines
! from shr_abort_abort. (This is for consistency with older code, from when these routines
! were defined in shr_sys_mod.)
public :: shr_sys_abort

! Fake sleep
public :: shr_sys_sleep

! Real flush
public :: shr_sys_flush

contains

subroutine shr_sys_sleep(sec)
  real(shr_kind_r8), intent(in) :: sec

  ! do nothing
end subroutine shr_sys_sleep

SUBROUTINE shr_sys_flush(unit)

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

   flush(unit)

END SUBROUTINE shr_sys_flush

end module shr_sys_mod
