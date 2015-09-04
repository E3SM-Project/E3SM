module clm_time_manager
  
  ! This is a stub for clm_time_manager. For a number of routines, it returns fixed values
  ! rather than going through all the time manager stuff (which would be awkward from a
  ! unit test).

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

contains

  integer function get_step_size()

    get_step_size = 1800

  end function get_step_size

  logical function is_first_step()

    is_first_step = .false.

  end function is_first_step

  integer function get_driver_start_ymd

    get_driver_start_ymd = 0

  end function get_driver_start_ymd

  integer function get_nstep

    get_nstep = 0

  end function get_nstep

  logical function is_end_curr_day

    is_end_curr_day = .false.

  end function is_end_curr_day

  subroutine get_start_date(yr, mon, day, tod)

    ! Return date components valid at beginning of initial run.

    ! Arguments
    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    yr = 0
    mon = 0
    day = 0
    tod = 0

  end subroutine get_start_date

  subroutine get_curr_date(yr, mon, day, tod, offset)
    
    !-----------------------------------------------------------------------------------------
    ! Return date components valid at end of current timestep with an optional
    ! offset (positive or negative) in seconds.
    !
    ! Return a fixed date, rather than going through the clm_time_manager, which would
    ! be awkward from a unit test

    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    yr = 1
    mon = 1
    day = 1
    tod = 0

  end subroutine get_curr_date

  function get_curr_yearfrac( offset )

    !---------------------------------------------------------------------------------
    ! Get the fractional position in the current year. This is 0 at midnight on Jan 1,
    ! and 1 at the end of Dec 31.
    !
    ! Return a fixed weight, rather than going through the clm_time_manager, which would
    ! be awkward from a unit test

    !
    ! Arguments
    real(r8) :: get_curr_yearfrac  ! function result
    
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds (ignored).

    real(r8), parameter :: fixed_weight = 0.75_r8

    get_curr_yearfrac = fixed_weight
  end function get_curr_yearfrac

end module clm_time_manager
