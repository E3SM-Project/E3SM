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
!  SVN:$Id: output.F90 46059 2013-04-16 22:02:01Z mlevy@ucar.edu $
!
! !USES:

   use kinds_mod
   use domain
!   use constants, only: 
!   use time_management, only: 
   use restart, only: write_restart, init_restart, lrestart_write
   use history, only: write_history, init_history
   use movie, only: write_movie, init_movie
   use overflows
   use overflow_type
   use tavg, only: write_tavg, init_tavg, final_tavg
   use timers, only: get_timer, timer_start, timer_stop

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: output_driver, &
             init_output,   &
             final_output


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
      timer_out,               &! timer for tavg
      timer_rest                ! timer for restart

   logical (log_kind), save :: &
      first_call = .true.       ! flag for initializing timers


!-----------------------------------------------------------------------
!
!  if this is the first call to output_driver, start some timers
!
!-----------------------------------------------------------------------


   if (first_call) then
      call get_timer(timer_out, 'OUTPUT'     ,nblocks_clinic,distrb_clinic%nprocs)
      call get_timer(timer_tavg,'OUTPUT TAVG',nblocks_clinic,distrb_clinic%nprocs)
      call get_timer(timer_rest,'OUTPUT REST',nblocks_clinic,distrb_clinic%nprocs)
      first_call = .false.
   endif

   call timer_start(timer_out)

!-----------------------------------------------------------------------
!
!  write history, movie files - the decision when to write
!  is internal to each routine  
!  write these first so that if I/O fails, no restart is written
!
!-----------------------------------------------------------------------

#if (defined _NOIO) 
! Insufficient memory to write history files on Blue Gene
#else
   call write_history
   call write_movie
#endif

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
!  restart must be written. 

!  note that lrestart_write is now a module variable, which allows
!  the overflows module to coordinate an overflows restart write
!
!-----------------------------------------------------------------------

   call timer_start(timer_tavg)
#if (defined _NOIO)
! Insufficient memory to restart tavg files on Blue Gene
#else
   call write_tavg(restart_type)
#endif
   call timer_stop (timer_tavg)

!-----------------------------------------------------------------------
!
!  write overflow restart file (coordinate with write_restart)
!
!-----------------------------------------------------------------------

    if ( lrestart_write .and. overflows_on .and. overflows_interactive ) then
       call ovf_write_restart
    endif

   call timer_stop(timer_out)
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
!BOP
! !IROUTINE: final_output
! !INTERFACE:

 subroutine final_output

! !DESCRIPTION:
!  Closes any files that are still open for output
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  Make sure tavg files are closed (note that history, movie, and restart
!  files are always closed after writing... tavg files are the only ones
!  that can be left open for multiple time steps)
!
!-----------------------------------------------------------------------

   call final_tavg

!-----------------------------------------------------------------------
!EOC

 end subroutine final_output


!***********************************************************************

 end module output

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
