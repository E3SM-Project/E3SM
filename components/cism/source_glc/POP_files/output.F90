!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module output

!BOP
! !MODULE: output
! !DESCRIPTION:
!  Contains necessary routines, variables for large model output
!  files - restart, history, movies, drifter, time average files.
!  This module is primarily a driver for the individual output
!  modules.
!
! !REVISION HISTORY:
!  SVN:$Id: output.F90 2290 2006-10-25 18:23:10Z njn01 $
!
! !USES:

   use kinds_mod
   use domain
!   use constants, only: 
!   use time_management, only: 
   use restart, only: write_restart, init_restart
   use history, only: write_history, init_history
   use movie, only: write_movie, init_movie
   use tavg, only: write_tavg, init_tavg
   use timers, only: get_timer, timer_start, timer_stop

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: output_driver, &
             init_output


!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: output_driver
! !INTERFACE:

   subroutine output_driver

! !DESCRIPTION:
!  This is the main driver routine for all large model output routines.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      restart_type          ! type of restart being written - used
                            ! to pass restart info to tavg routines

   integer (int_kind), save :: &
      timer_tavg,              &! timer for tavg
      timer_rest                ! timer for restart

   logical (log_kind), save :: &
      first_call = .true.       ! flag for initializing timers


!-----------------------------------------------------------------------
!
!  if this is the first call to output_driver, start some timers
!
!-----------------------------------------------------------------------

   if (first_call) then
      call get_timer(timer_tavg,'OUTPUT TAVG',nblocks_clinic,distrb_clinic%nprocs)
      call get_timer(timer_rest,'OUTPUT REST',nblocks_clinic,distrb_clinic%nprocs)
      first_call = .false.
   endif

!-----------------------------------------------------------------------
!
!  write history, movie files - the decision when to write
!  is internal to each routine  
!  write these first so that if I/O fails, no restart is written
!
!-----------------------------------------------------------------------

   call write_history
   call write_movie

!-----------------------------------------------------------------------
!
!  check for restart and write restart if required
!
!-----------------------------------------------------------------------

   call timer_start(timer_rest)
   call write_restart(restart_type)
   call timer_stop (timer_rest)

!-----------------------------------------------------------------------
!
!  write tavg - the decision when to write
!  is internal to routine except for notifying tavg that a 
!  restart must be written
!
!-----------------------------------------------------------------------

   call timer_start(timer_tavg)
   call write_tavg(restart_type)
   call timer_stop (timer_tavg)

!-----------------------------------------------------------------------
!EOC

 end subroutine output_driver

!***********************************************************************
!BOP
! !IROUTINE: init_output
! !INTERFACE:

 subroutine init_output

! !DESCRIPTION:
!  Initializes frequency of output and filenames for
!  various files by calling individual initialization routines
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  call individual init routines
!
!-----------------------------------------------------------------------

   call init_restart
   call init_history
   call init_movie
   call init_tavg

!-----------------------------------------------------------------------
!EOC

 end subroutine init_output


!***********************************************************************

 end module output

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
