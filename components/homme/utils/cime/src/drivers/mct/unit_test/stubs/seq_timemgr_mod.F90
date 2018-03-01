module seq_timemgr_mod

  ! Stub for routines from seq_timemgr_mod that are needed by other modules built by the
  ! unit tests.

  implicit none
  private

  public :: seq_timemgr_pause_active

contains

  logical function seq_timemgr_pause_active()
    ! Stub for seq_timemgr_pause_active - always returns .false.

    seq_timemgr_pause_active = .false.
  end function seq_timemgr_pause_active

end module seq_timemgr_mod
