!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_InitMod

!BOP
! !MODULE: POP_InitMod
! !DESCRIPTION:
!  This module contains the POP initialization method and initializes
!  everything needed by a POP simulation.  Primarily it is a driver
!  that calls individual initialization routines for each POP module.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_InitMod.F90 808 2006-04-28 17:06:38Z njn01 $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use kinds_mod, only: int_kind
   use initial
   use domain, only: distrb_clinic
   use timers, only: get_timer
   use time_management, only: access_time_flag

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_Initialize1, POP_Initialize2

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------
   integer (int_kind), public :: &
      stop_now,                  &! flag id for stop_now flag
      cpl_ts,                    &! flag id for coupled timestep flag
      timer_total,               &! timer for entire run phase
      nscan


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_Initialize1
! !INTERFACE:

 subroutine POP_Initialize1(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  initialize return flag
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!  call pop initialization routines
!
!-----------------------------------------------------------------------

   call pop_init_phase1(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_Initialize1: error in pop_init_phase1')
      return
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize1

!***********************************************************************

!BOP
! !IROUTINE: POP_Initialize2
! !INTERFACE:

 subroutine POP_Initialize2(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize return flag
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!  complete pop initialization process
!
!-----------------------------------------------------------------------

   call pop_init_phase2(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_Initialize1: error in pop_init_phase2')
      return
   endif

!-----------------------------------------------------------------------
!
!  initialize variables used by the pop driver
!
!-----------------------------------------------------------------------
   nscan = 0

!-----------------------------------------------------------------------
!
!  initialize driver-level flags and timers
!
!-----------------------------------------------------------------------
   call access_time_flag ('stop_now', stop_now)
   call access_time_flag ('coupled_ts', cpl_ts)

   call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize2

!***********************************************************************

 end module POP_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
