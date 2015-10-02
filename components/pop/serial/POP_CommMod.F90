!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module POP_CommMod

! !MODULE: POP_CommMod
! !DESCRIPTION:
!  This module contains necessary routines and variables to support
!  other parallel communication modules in POP.  In particular, this
!  module contains communicators, tags, task ids and other necessary
!  information and the routines to set them up.  In addition, several
!  utility routines for setting up the communication environment are
!  included.  For this serial version, most of these routines simply
!  set dummy values and perform no operations.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_CommMod.F90 52 2007-02-26 20:08:12Z pwjones $
!  2006-07-10: Phil Jones
!              add new communication module with new naming convention
!              also contains new coupler code based on merge between
!                 POP 2.0 code and NCAR CCSM POP from Nancy Norton
!
! !USES:

   use POP_KindsMod

#ifdef CCSMCOUPLED
   use cpl_interface_mod, only : cpl_interface_init,cpl_interface_finalize
   use cpl_fields_mod, only : cpl_fields_ocnname
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: POP_CommInit,                    &
              POP_CommInitMessageEnvironment,  &
              POP_CommExitMessageEnvironment,  &
              POP_CommAbortMessageEnvironment, &
              POP_CommGetNumProcs,             &
              POP_CommCreateCommunicator,      &
              POP_Barrier

! !PUBLIC DATA MEMBERS:

   integer (POP_i4), public :: &
      POP_communicator,         &! MPI communicator for ocn comms
      POP_mpiR8,                &! MPI type for r8
      POP_mpiR4,                &! MPI type for r4
      POP_myTask,               &! MPI task number for this task
      POP_masterTask             ! MPI task number for master task

   integer (POP_i4), parameter, public :: &
      POP_mpitagBndy2d         = 1,       &! MPI tags for various
      POP_mpitagGS             = 1000      ! communication patterns

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_CommInit
! !INTERFACE:

 subroutine POP_CommInit

! !DESCRIPTION:
!  This routine sets up communication environment and defines ocean
!  communicator.
!
! !REVISION HISTORY:
!  same as module
!

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  create communicator for internal ocean communications
!  this process varies for different coupling styles and is often
!  determined by a library call.  These are determined by a preprocessor
!  directive.  It is assumed that the MPI initialization occurs within
!  these calls or that POP_CommInitMessageEnvironment is called from
!  the driver routine.
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  call CCSM coupler routine to return communicator
!
!-----------------------------------------------------------------------


   call cpl_interface_init(cpl_fields_ocnname, POP_communicator)

#else
!-----------------------------------------------------------------------
!
!  when not coupled, simply set a dummy value
!
!-----------------------------------------------------------------------

   POP_communicator = 0

#endif

!-----------------------------------------------------------------------
!
!  initialize variables with dummy value for the serial case
!
!-----------------------------------------------------------------------

   POP_masterTask = 0
   POP_myTask     = 0
   POP_mpiR8      = 0
   POP_mpiR4      = 0

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommInit

!***********************************************************************
!BOP
! !IROUTINE: POP_CommGetNumProcs
! !INTERFACE:

 function POP_CommGetNumProcs(communicator)

! !DESCRIPTION:
!  This function returns the number of processor assigned to
!  a given communicator.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      communicator         ! communicator to query for num processors

! !OUTPUT PARAMETERS:

   integer (POP_i4) :: &
      POP_CommGetNumProcs  ! number of processors in communicator

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  always return one for serial case
!
!-----------------------------------------------------------------------

   POP_CommGetNumProcs = 1

!-----------------------------------------------------------------------
!EOC

 end function POP_CommGetNumProcs

!***********************************************************************
!BOP
! !IROUTINE: POP_CommInitMessageEnvironment
! !INTERFACE:

 subroutine POP_CommInitMessageEnvironment

! !DESCRIPTION:
!  This routine initializes the message environment.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  this routine does nothing in serial case
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommInitMessageEnvironment

!***********************************************************************
!BOP
! !IROUTINE: POP_CommExitMessageEnvironment
! !INTERFACE:

 subroutine POP_CommExitMessageEnvironment

! !DESCRIPTION:
!  This routine exits the message environment properly when model
!  stops.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for coupled model, call cpl routine
!  otherwise, serial case does nothing
!
!-----------------------------------------------------------------------


#ifdef CCSMCOUPLED
   call cpl_interface_finalize(cpl_fields_ocnname)
#else

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommExitMessageEnvironment

!***********************************************************************
!BOP
! !IROUTINE: POP_CommAbortMessageEnvironment
! !INTERFACE:

 subroutine POP_CommAbortMessageEnvironment

! !DESCRIPTION:
!  This routine aborts the message environment when model stops.
!  It will attempt to abort the entire MPI COMM WORLD.
!
! !REVISION HISTORY:
!  same as module

#ifdef CCSMCOUPLED
! !INCLUDES

   include 'mpif.h'
#endif

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  !MPI error flag

!-----------------------------------------------------------------------
!
!  call cpl routines to abort if coupled
!  otherwise, serial case does nothing
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   ierr = 13
   call MPI_ABORT(0,ierr)
   call cpl_interface_finalize(cpl_fields_ocnname)
#else
   ierr = 0
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommAbortMessageEnvironment

!***********************************************************************
!BOP
! !IROUTINE: POP_CommCreateCommunicator
! !INTERFACE:

 subroutine POP_CommCreateCommunicator(newCommunicator, numProcs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ocean communicator.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      numProcs          ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      newCommunicator   ! new communicator for this distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set dummy value for serial case
!
!-----------------------------------------------------------------------

   newCommunicator = 0

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommCreateCommunicator

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: POP_Barrier
! !INTERFACE:

 subroutine POP_Barrier

! !DESCRIPTION:
!   Because this is serial code this routine does nothing
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Barrier


 end module POP_CommMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
