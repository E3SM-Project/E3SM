module clm_time_manager
  
  ! This is a mock replacement for clm_time_manager. It mocks out get_curr_yearfrac,
  ! returning a fixed weight rather than going through all the time manager stuff (which
  ! would be awkward from a unit test).

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

contains

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
