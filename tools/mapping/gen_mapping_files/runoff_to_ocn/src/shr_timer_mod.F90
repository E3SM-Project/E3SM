module shr_timer_mod

   !----------------------------------------------------------------------------
   !
   ! routines that support multiple CPU timers via F90 intrinsics
   !
   ! Note:
   ! o if   an operation is requested on an invalid timer number n
   !   then nothing is done in a routine
   ! o if   more than max_timers are requested,
   !   then timer n=max_timers is "overloaded" and becomes invalid/undefined
   !
   ! * cpp if-defs were introduced in 2005 to work-around a bug in the ORNL Cray
   !   X1 F90 intrinsic system_clock() function -- ideally this Cray bug would be
   !   fixed and cpp if-defs would be unnecessary and removed.
   !
   ! !REVISION HISTORY:
   !    2005-??-?? - added workaround for Cray F90 bug, mods by Cray/ORNL
   !    2000-??-?? - 1st version by B. Kauffman
   !----------------------------------------------------------------------------

   use shr_kind_mod

   implicit none

   private  ! restricted access

   public  :: shr_timer_init , shr_timer_get
   public  :: shr_timer_start, shr_timer_stop
   public  :: shr_timer_print, shr_timer_print_all
   public  :: shr_timer_check, shr_timer_check_all
   public  :: shr_timer_zero , shr_timer_zero_all
   public  :: shr_timer_free , shr_timer_free_all
   public  :: shr_timer_sleep

   integer(SHR_KIND_IN),parameter :: stat_free    = 0  ! timer status constants
   integer(SHR_KIND_IN),parameter :: stat_inuse   = 1
   integer(SHR_KIND_IN),parameter :: stat_started = 2
   integer(SHR_KIND_IN),parameter :: stat_stopped = 3
   integer(SHR_KIND_IN),parameter :: max_timers   = 200 ! max number of timers

   integer(SHR_KIND_IN) :: status (max_timers) ! status of each timer
   !----------------------------------------------------------------------------
   ! the following ifdef circumvents a bug in the X1 system_clock function
   !----------------------------------------------------------------------------
#if (defined UNICOSMP)
   integer(kind=8)      :: cycles1(max_timers) ! cycle number at timer start
   integer(kind=8)      :: cycles2(max_timers) ! cycle number at timer stop
#else
   integer(SHR_KIND_IN) :: cycles1(max_timers) ! cycle number at timer start
   integer(SHR_KIND_IN) :: cycles2(max_timers) ! cycle number at timer stop
#endif
   integer(SHR_KIND_IN) :: cycles_max = -1     ! max cycles before wrapping
   character   (len=80) :: name   (max_timers) ! name assigned to each timer
   real   (SHR_KIND_R8) :: dt     (max_timers) ! accumulated time
   integer(SHR_KIND_IN) :: calls  (max_timers) ! # of samples in accumulation
   real   (SHR_KIND_R8) :: clock_rate          ! clock_rate: seconds per cycle

   save

!===============================================================================
   contains
!===============================================================================

subroutine shr_timer_init

   !----- local -----
   integer(SHR_KIND_IN) :: cycles ! count rate return by system clock
#if (defined UNICOSMP)
   integer(kind=8) :: irtc_rate
#endif

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_init) ',a,i5)"

!-------------------------------------------------------------------------------
! This routine initializes:
! 1) values in all timer array locations
! 2) machine parameters necessary for computing cpu time from F90 intrinsics.
!    F90 intrinsic: system_clock(count_rate=cycles, count_max=cycles_max)
!-------------------------------------------------------------------------------

   call shr_timer_free_all

#if (defined UNICOSMP)
   cycles = irtc_rate()
#else
   call system_clock(count_rate=cycles, count_max=cycles_max)
#endif

   if (cycles /= 0) then
     clock_rate = 1.0_SHR_KIND_R8/real(cycles,SHR_KIND_R8)
   else
     clock_rate = 0._SHR_KIND_R8
     write(6,F00) 'ERROR: no system clock available'
   endif

end subroutine shr_timer_init

!===============================================================================

subroutine shr_timer_get(n, str)

   !----- arguments -----
   integer(SHR_KIND_IN),intent(out) :: n    ! timer number
   character (*)       ,intent( in) :: str  ! text string with timer name

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_get) ',a,i5)"

!-----------------------------------------------------------------------
!  search for next free timer
!-----------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_free) then
       status(n) = stat_inuse
       name  (n) = str
       calls (n) = 0
       return
     endif
   end do

   n=max_timers
   name  (n) = "<invalid - undefined - overloaded>"
   write(6,F00) 'ERROR: exceeded maximum number of timers'

end subroutine shr_timer_get

!===============================================================================

subroutine shr_timer_start(n)

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in) :: n      ! timer number

   !----- local -----
#if (defined UNICOSMP)
   integer(kind=8) :: irtc
#endif

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_start) ',a,i5)"

!-----------------------------------------------------------------------
!  This routine starts a given timer.
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) call shr_timer_stop(n)

     status(n) = stat_started
#if (defined UNICOSMP)
     cycles1(n) = irtc()
#else
     call system_clock(count=cycles1(n))
#endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_start

!===============================================================================

subroutine shr_timer_stop(n)

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in) :: n  ! timer number

   !----- local -----
   real (SHR_KIND_R8) :: elapse      ! elapsed time returned by system counter
#if (defined UNICOSMP)
   integer(kind=8) :: irtc
#endif

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_stop) ',a,i5)"

!-------------------------------------------------------------------------------
!  This routine stops a given timer, checks for cycle wrapping, computes the
!  elapsed time, and accumulates the elapsed time in the dt(n) array
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if ( status(n) == stat_started) then
#if (defined UNICOSMP)
       cycles2(n) = irtc()
#else
       call system_clock(count=cycles2(n))
#endif
       if (cycles2(n) >= cycles1(n)) then
         dt(n) = dt(n) + clock_rate*(cycles2(n) - cycles1(n))
       else
         dt(n) = dt(n) + clock_rate*(cycles_max + cycles2(n) - cycles1(n))
       endif
       calls (n) = calls(n) + 1
       status(n) = stat_stopped
     end if
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_stop

!===============================================================================

subroutine shr_timer_print(n)

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in) :: n     ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_print) ',a,i5)"
   character(len=*),parameter :: F01 = "('(shr_timer_print) timer',i3,&
   &                                     ':',i8,' calls,',f10.3,'s, id: ',a)"
!-------------------------------------------------------------------------------
!  prints the accumulated time for a given timer
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call shr_timer_stop(n)
       write (6,F01) n,calls(n),dt(n),trim(name(n))
       call shr_timer_start(n)
     else
       write (6,F01) n,calls(n),dt(n),trim(name(n))
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_print

!===============================================================================

subroutine shr_timer_print_all

   !----- local -----
   integer(SHR_KIND_IN) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_print_all) ',a,i5)"

!-------------------------------------------------------------------------------
!  prints accumulated time for all timers in use
!-------------------------------------------------------------------------------

   write(6,F00) 'print all timing info:'

   do n=1,max_timers
     if (status(n) /= stat_free) call shr_timer_print(n)
   end do

end subroutine shr_timer_print_all

!===============================================================================

subroutine shr_timer_zero(n)

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in) :: n       ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_zero) ',a,i5)"

!-------------------------------------------------------------------------------
!  This routine resets a given timer.
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     dt(n) = 0.0_SHR_KIND_R8
     calls(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_zero

!===============================================================================

subroutine shr_timer_zero_all

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_zero_all) ',a,i5)"

!-------------------------------------------------------------------------------
!  This routine resets all timers.
!-------------------------------------------------------------------------------

   dt = 0.0_SHR_KIND_R8
   calls = 0

end subroutine shr_timer_zero_all

!===============================================================================

subroutine shr_timer_check(n)

   !----- arguments -----
   integer(SHR_KIND_IN), intent(in) ::  n   ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_check) ',a,i5)"

!-------------------------------------------------------------------------------
!  This routine checks a given timer.  This is primarily used to
!  periodically accumulate time in the timer to prevent timer cycles
!  from wrapping around max_cycles.
!-------------------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     if (status(n) == stat_started) then
       call shr_timer_stop (n)
       call shr_timer_start(n)
     endif
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_check

!===============================================================================

subroutine shr_timer_check_all

   !----- local -----
   integer(SHR_KIND_IN) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_check_all) ',a,i5)"

!-------------------------------------------------------------------------------
!  Call shr_timer_check for all timers in use
!-------------------------------------------------------------------------------

   do n=1,max_timers
     if (status(n) == stat_started) then
       call shr_timer_stop (n)
       call shr_timer_start(n)
     endif
   end do

end subroutine shr_timer_check_all

!===============================================================================

subroutine shr_timer_free(n)

   !----- arguments -----
   integer(SHR_KIND_IN),intent(in) :: n    ! timer number

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_free) ',a,i5)"

!-----------------------------------------------------------------------
!  initialize/free all timer array values
!-----------------------------------------------------------------------

   if ( n>0 .and. n<=max_timers) then
     status (n) = stat_free
     name   (n) = "<invalid - undefined>"
     dt     (n) = 0.0_SHR_KIND_R8
     cycles1(n) = 0
     cycles2(n) = 0
   else
     write(6,F00) 'ERROR: invalid timer number: ',n
   end if

end subroutine shr_timer_free

!===============================================================================

subroutine shr_timer_free_all

   !----- local -----
   integer(SHR_KIND_IN) :: n

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(shr_timer_free_all) ',a,i5)"

!-------------------------------------------------------------------------------
!  initialize/free all timer array values
!-------------------------------------------------------------------------------

   do n=1,max_timers
     call shr_timer_free(n)
   end do

end subroutine shr_timer_free_all

!===============================================================================

subroutine shr_timer_sleep(sec)

   use shr_sys_mod     ! share system calls (namely, shr_sys_sleep)

   !----- local -----
   real   (SHR_KIND_R8),intent(in) :: sec  ! number of seconds to sleep

!-------------------------------------------------------------------------------
! Sleep for approximately sec seconds
!
! Note: sleep is typically a system call, hence it is implemented in
!       shr_sys_mod, although it probably would only be used in a timing
!       context, which is why there is a shr_timer_* wrapper provided here.
!-------------------------------------------------------------------------------

   call shr_sys_sleep(sec)

end subroutine shr_timer_sleep

!===============================================================================
end module shr_timer_mod
!===============================================================================
