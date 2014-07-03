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
!  included.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_CommMod.F90 72 2007-10-05 21:01:46Z pwjones $
!  2006-07-10: Phil Jones
!              add new communication module with new naming conventions
!              also adds new CCSM coupler interface based on merge
!                 between POP 2.0 code and CCSM POP from Nancy Norton
!
! !USES:

   use POP_KindsMod

#ifdef CCSMCOUPLED
   use ocn_communicator, only: mpi_communicator_ocn
#endif

   implicit none
   private
   save

   include 'mpif.h'

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
      POP_mpiR16,               &! MPI type for r16
      POP_mpiR8,                &! MPI type for r8
      POP_mpiR4,                &! MPI type for r4
      POP_myTask,               &! MPI task number for this task
      POP_masterTask             ! MPI task number for master task

   integer (POP_i4), parameter, public :: &
      POP_mpitagHalo           = 1,       &! MPI tags for various
      POP_mpitagRedist         = 1000      ! communication patterns

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
!  This routine sets up MPI environment and defines ocean
!  communicator.
!
! !REVISION HISTORY:
!  same as module
!
! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! MPI error flag

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
   POP_communicator = mpi_communicator_ocn

#else
!-----------------------------------------------------------------------
!
!  when not coupled, simply duplicate global communicator
!
!-----------------------------------------------------------------------

   call MPI_COMM_DUP(MPI_COMM_WORLD, POP_communicator, ierr)

#endif

!-----------------------------------------------------------------------
!
!  determine task ids
!
!-----------------------------------------------------------------------

   POP_masterTask = 0
   call MPI_COMM_RANK  (POP_communicator, POP_myTask, ierr)

!-----------------------------------------------------------------------
!
!  On some machines the MPI implementation makes some assumptions about
!  these data types, so these are chosen to try and choose the
!  appropriate kind.
!
!-----------------------------------------------------------------------

   POP_mpiR16 = MPI_REAL16
   POP_mpiR8  = MPI_REAL8
   POP_mpiR4  = MPI_REAL4

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
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr

!-----------------------------------------------------------------------

   call MPI_COMM_SIZE(communicator, POP_CommGetNumProcs, ierr)

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

! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! MPI error flag

!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   ! initialized by cpl library routines
#else
   call MPI_INIT(ierr)
#endif

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

! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! MPI error flag

!-----------------------------------------------------------------------

#ifndef CCSMCOUPLED
   call MPI_FINALIZE(ierr)
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

! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: errorCode, ierr  !MPI error flag

!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
!  call MPI_BARRIER(POP_Communicator,ierr)
   ierr = 13
   call MPI_ABORT(MPI_COMM_WORLD,errorCode, ierr)
#else
   call MPI_BARRIER(POP_Communicator, ierr)
   call MPI_ABORT(MPI_COMM_WORLD, errorCode, ierr)
   call MPI_FINALIZE(ierr)
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

! !INCLUDES:

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
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     MPI_GROUP_OCN,         &! group of processors assigned to ocn
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (POP_i4) :: &
     ierr                    ! error flag for MPI comms

   integer (POP_i4), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP (POP_Communicator, MPI_GROUP_OCN, ierr)

   range(1) = 0
   range(2) = numProcs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in POP_Communicator
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_OCN, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (POP_Communicator, MPI_GROUP_NEW,  &
                         newCommunicator, ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommCreateCommunicator

!***********************************************************************
!BOP
! !IROUTINE: POP_Barrier
! !INTERFACE:

 subroutine POP_Barrier

! !DESCRIPTION:
!  This routine performs a barrier.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer(POP_i4) :: ierr

    call MPI_Barrier(POP_Communicator,ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Barrier


 end module POP_CommMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
