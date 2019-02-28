!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module uses F90 cpu time routines to allowing setting of
!     multiple CPU timers.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: timers.f,v 1.2 2000/04/19 21:56:26 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!     This code has been modified from the version available from 
!     Los Alamos National Laboratory, for the purpose of running it
!     within WW3.
!
!***********************************************************************

      module scrip_timers

!-----------------------------------------------------------------------

      use SCRIP_KindsMod ! defines common data types

      implicit none

      integer (SCRIP_i4), parameter ::  
     &     max_timers = 99  ! max number of timers allowed

      integer (SCRIP_i4), save :: 
     &     cycles_max       ! max value of clock allowed by system

      integer (SCRIP_i4), dimension(max_timers), save :: 
     &     cycles1,         ! cycle number at start for each timer
     &     cycles2          ! cycle number at stop  for each timer

      real (SCRIP_r4), save ::  
     &     clock_rate       ! clock_rate in seconds for each cycle

      real (SCRIP_r4), dimension(max_timers), save ::  
     &     cputime          ! accumulated cpu time in each timer

      character (len=8), dimension(max_timers), save ::  
     &     status           ! timer status string

!***********************************************************************

      contains

!***********************************************************************

      subroutine timer_check(timer)

!-----------------------------------------------------------------------
!
!     This routine checks a given timer.  This is primarily used to
!     periodically accumulate time in the timer to prevent timer cycles
!     from wrapping around max_cycles.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------

      if (status(timer) .eq. 'running') then
        call timer_stop (timer)
        call timer_start(timer)
      endif

!-----------------------------------------------------------------------

      end subroutine timer_check

!***********************************************************************

      subroutine timer_clear(timer)

!-----------------------------------------------------------------------
!
!     This routine resets a given timer.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------

      cputime(timer) = 0.0_SCRIP_r4  ! clear the timer

!-----------------------------------------------------------------------

      end subroutine timer_clear

!***********************************************************************

      function timer_get(timer)

!-----------------------------------------------------------------------
!
!     This routine returns the result of a given timer.  This can be
!     called instead of timer_print so that the calling routine can 
!     print it in desired format.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------
!
!     Output Variables:
!
!-----------------------------------------------------------------------

      real (SCRIP_r4) ::  
     &     timer_get   ! accumulated cputime in given timer

!-----------------------------------------------------------------------

      if (status(timer) .eq. 'stopped') then
        timer_get = cputime(timer)
      else
        call timer_stop(timer)
        timer_get = cputime(timer)
        call timer_start(timer)
      endif

!-----------------------------------------------------------------------

      end function timer_get

!***********************************************************************

      subroutine timer_print(timer)

!-----------------------------------------------------------------------
!
!     This routine prints the accumulated cpu time in given timer.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------

      !--- 
      !--- print the cputime accumulated for timer 
      !--- make sure timer is stopped
      !---

      if (status(timer) .eq. 'stopped') then
        write(*,"(' CPU time for timer',i3,':',1p,e16.8)")  
     &       timer,cputime(timer)
      else
        call timer_stop(timer)
        write(*,"(' CPU time for timer',i3,':',1p,e16.8)")  
     &       timer,cputime(timer)
        call timer_start(timer)
      endif

!-----------------------------------------------------------------------

      end subroutine timer_print

!***********************************************************************

      subroutine timer_start(timer)

!-----------------------------------------------------------------------
!
!     This routine starts a given timer.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------

      !---
      !--- Start the timer and change timer status.
      !---

      if (status(timer) .eq. 'stopped') then
        call system_clock(count=cycles1(timer))
        status(timer) = 'running'
      endif

!-----------------------------------------------------------------------

      end subroutine timer_start

!***********************************************************************

      subroutine timer_stop(timer)

!-----------------------------------------------------------------------
!
!     This routine stops a given timer.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input Variables:
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::  
     &    timer            ! timer number

!-----------------------------------------------------------------------

      if (status(timer) .eq. 'running') then

        !---
        !--- Stop the desired timer.
        !---

        call system_clock(count=cycles2(timer))

        !---
        !--- check and correct for cycle wrapping
        !---

        if (cycles2(timer) .ge. cycles1(timer)) then
          cputime(timer) = cputime(timer) + clock_rate*  
     &                     (cycles2(timer) - cycles1(timer))
        else
          cputime(timer) = cputime(timer) + clock_rate*  
     &                (cycles2(timer) - cycles1(timer) + cycles_max)
        endif

        !---
        !--- Change timer status.
        !---

        status(timer)='stopped'

      endif

!-----------------------------------------------------------------------

      end subroutine timer_stop

!***********************************************************************

      subroutine timers_init

!-----------------------------------------------------------------------
!
!     This routine initializes some machine parameters necessary for
!     computing cpu time from F90 intrinsics.
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: cycles ! count rate return by sys_clock

!-----------------------------------------------------------------------

      !---
      !--- Initialize timer arrays and clock_rate.
      !---

      clock_rate = 0.0_SCRIP_r4
      cycles1    = 0
      cycles2    = 0
      cputime    = 0.0_SCRIP_r4
      status     = 'stopped'

      !---
      !--- Call F90 intrinsic system_clock to determine clock rate
      !--- and maximum cycles.  If no clock available, print message.
      !---

      call system_clock(count_rate=cycles, count_max=cycles_max)

      if (cycles /= 0) then
        clock_rate = 1.0_SCRIP_r4/real(cycles)
      else
        clock_rate = 0.0_SCRIP_r4
        print *, '--- No system clock available ---'
      endif

!-----------------------------------------------------------------------

      end subroutine timers_init

!***********************************************************************

      end module scrip_timers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
