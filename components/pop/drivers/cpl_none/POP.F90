!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !ROUTINE: POP
! !INTERFACE:

 program POP

! !DESCRIPTION:
!  This is the main driver for the Parallel Ocean Program (POP).
!
! !REVISION HISTORY:
!  CVS:$Id: POP.F90,v 1.8 2003/01/28 23:21:19 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_InitMod

   use kinds_mod, only: int_kind, r8
   use communicate, only: my_task, master_task
!   use constants, only: 
   use domain, only: distrb_clinic
   use timers, only: timer_print_all, get_timer, timer_start, timer_stop
   use time_management, only: get_time_flag_id, check_time_flag,    &
       nsteps_run, stdout, exit_pop, override_time_flag
   use step_mod, only: step
   use diagnostics, only: check_KE
   use output, only: output_driver
   use exit_mod, only: sigAbort, sigExit
!   use io, only: 

   implicit none

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode          ! error flag

   integer (int_kind) :: &
      timer_step,        &! timer number for step
      fstop_now           ! flag id for stop_now flag

!-----------------------------------------------------------------------
!
!  initialize the model run
!
!-----------------------------------------------------------------------

   call POP_CommInitMessageEnvironment

   errorCode = POP_Success

   call POP_Initialize(errorCode)

   fstop_now = get_time_flag_id('stop_now')

!-----------------------------------------------------------------------
!
!  start up the main timer
!
!-----------------------------------------------------------------------

   call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)
   call timer_start(timer_total)

   call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!
!  advance the model in time
!
!-----------------------------------------------------------------------

   do while (.not. check_time_flag(fstop_now) .and. &
             errorCode == POP_Success)

      call timer_start(timer_step)
      call step(errorCode)
      call timer_stop(timer_step)

      !***
      !*** exit if energy is blowing
      !***

      if (check_KE(100.0_r8)) then
         call override_time_flag(fstop_now,value=.true.)
         call output_driver
         call exit_POP(sigAbort,'ERROR: k.e. > 100 ')
      endif

!-----------------------------------------------------------------------
!
!     write restart dumps and archiving
!
!-----------------------------------------------------------------------

      call output_driver

   enddo

!-----------------------------------------------------------------------
!
!  print timing information and clean up various environments if 
!  they have been used
!
!-----------------------------------------------------------------------

   call timer_stop(timer_total)
   call timer_print_all(stats=.true.)

   call POP_ErrorPrint(errorCode, printTask=POP_masterTask)

   call exit_POP(sigExit,'Successful completion of POP run')

!-----------------------------------------------------------------------
!EOC

 end program POP

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
