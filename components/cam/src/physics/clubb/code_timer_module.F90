! $Id: code_timer_module.F90 8106 2016-05-17 23:29:01Z raut@uwm.edu $
module code_timer_module

! Description:
!   This module contains a diagnostic timer utility that can be used
!   to time a piece of code.

  implicit none

  private ! Set default scope

  ! A timer!!
  type timer_t
    real :: time_elapsed        ! Time elapsed [sec]
    real :: secstart            ! Timer starting time
  end type timer_t

  public :: timer_t, timer_start, timer_stop

  contains

  !-----------------------------------------------------------------------
  subroutine timer_start( timer )

  ! Description:
  !   Starts the timer

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variables
    type(timer_t), intent(inout) :: timer

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    call cpu_time( timer%secstart )
    return
  end subroutine timer_start
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_stop( timer )

  ! Description:
  !   Stops the timer

  ! References:
  !   None
  !-----------------------------------------------------------------------
    implicit none

    ! Input/Output Variables
    type(timer_t), intent(inout) :: timer

    ! Local Variables
    real :: secend

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    call cpu_time( secend )


    timer%time_elapsed = timer%time_elapsed + (secend - timer%secstart)
    timer%secstart = 0.0

    return
  end subroutine timer_stop
  !-----------------------------------------------------------------------

end module code_timer_module
