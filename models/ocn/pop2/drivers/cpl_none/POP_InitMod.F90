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
!  SVN:$Id: POP_InitMod.F90 8528 2008-01-15 01:49:19Z dennis $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use kinds_mod, only: int_kind
   use initial
   use domain, only: distrb_clinic
   use timers, only: get_timer
   use time_management, only: get_time_flag_id

#ifdef coupled
   use POP_CouplingMod, only: pop_init_coupler_comm, pop_send_to_coupler, irbuf
   use cpl_interface_mod, only: cpl_interface_ibufRecv
   use cpl_fields_mod, only: cpl_fields_cplname
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_Initialize, POP_Initialize1, POP_Initialize2,  &
             POP_Initialize_coupling
             

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
! !IROUTINE: POP_Initialize 
! !INTERFACE:

 subroutine POP_Initialize(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines, with coupling
!  calls inbetween initialization calls. When invoking CCSM with coupling
!  at the top, the routines called in this subroutine will be accessed 
!  elsewhere, and this routine will not be called.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
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
!  call pop initialization routines in two stages, with coupling inbetween
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_Initialize1(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_Initialize: error in POP_Initialize1')
      return
   endif

!-----------------------------------------------------------------------
!
!  exchange initial information with coupler
!
!-----------------------------------------------------------------------

   call POP_Initialize_coupling

!-----------------------------------------------------------------------
!
!  complete the initialiation of the model
!
!-----------------------------------------------------------------------

   call POP_Initialize2(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_Initialize: error in POP_Initialize2')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize 

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

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
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

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
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
   stop_now  = get_time_flag_id('stop_now')
   cpl_ts    = get_time_flag_id('coupled_ts')

   call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize2

!***********************************************************************

!BOP
! !IROUTINE: POP_Initialize_coupling
! !INTERFACE:

 subroutine POP_Initialize_coupling

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
! !USES


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   
   logical (POP_Logical),save ::  &
      coupled_ts = .false.      ! coupled_ts is false at initialization


#if coupled

!-----------------------------------------------------------------------
!
!  exchange initial information with coupler
!
!-----------------------------------------------------------------------

   call pop_init_coupler_comm 

!-----------------------------------------------------------------------
!
!  receive initial message from coupler
!
!-----------------------------------------------------------------------

   call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

!-----------------------------------------------------------------------
!
!  send initial state information to the coupler
!
!-----------------------------------------------------------------------

   call pop_send_to_coupler(coupled_ts)

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize_coupling

 end module POP_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
