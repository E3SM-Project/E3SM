!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module communicate

! !MODULE: communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.
!
! !REVISION HISTORY:
!  SVN:$Id: communicate.F90 12674 2008-10-31 22:21:32Z njn01 $
!  2006-07-10: Phil Jones
!              edited to use new POP comm module - this module
!                 now for back compatibility only 
!
! !USES:

   use kinds_mod
   use POP_CommMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_communicate,          &
              exit_message_environment,  &
              abort_message_environment, &
              get_num_procs,             &
              create_communicator

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      MPI_COMM_OCN,             &! MPI communicator for ocn comms
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

   integer (int_kind), parameter, public :: &
      mpitag_bndy_2d        = 1,    &! MPI tags for various
      mpitag_bndy_3d        = 2,    &! communication patterns
      mpitag_gs             = 1000   ! communication patterns

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_communicate
! !INTERFACE:

 subroutine init_communicate

! !DESCRIPTION:
!  This routine sets up MPI environment and defines ocean
!  communicator.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  initialize communication routines and variables  
!  this interface for back compatibility only and assume POP_CommInit
!  and POP_CommInitMessageEnvironment have already been called
!
!-----------------------------------------------------------------------

   call POP_CommInit

   MPI_COMM_OCN = POP_communicator
   master_task  = POP_masterTask
   my_task      = POP_myTask
   MPI_DBL      = POP_mpiR8

!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: get_num_procs
! !INTERFACE:

 function get_num_procs()

! !DESCRIPTION:
!  This function returns the number of processor assigned to
!  MPI_COMM_OCN
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------

   get_num_procs = POP_CommGetNumProcs(POP_communicator)

!-----------------------------------------------------------------------
!EOC

 end function get_num_procs

!***********************************************************************
!BOP
! !IROUTINE: exit_message_environment
! !INTERFACE:

 subroutine exit_message_environment(ierr)

! !DESCRIPTION:
!  This routine exits the message environment properly when model
!  stops.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'   ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr   ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

   ierr = 0
   call POP_CommExitMessageEnvironment
 
!-----------------------------------------------------------------------
!EOC

 end subroutine exit_message_environment

!***********************************************************************
!BOP
! !IROUTINE: abort_message_environment
! !INTERFACE:

 subroutine abort_message_environment(ierr)

! !DESCRIPTION:
!  This routine aborts the message environment when model stops.
!  It will attempt to abort the entire MPI COMM WORLD.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'   ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr   ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

   ierr = 0
   call POP_CommAbortMessageEnvironment
   
!-----------------------------------------------------------------------
!EOC

 end subroutine abort_message_environment

!***********************************************************************
!BOP
! !IROUTINE: create_communicator
! !INTERFACE:

 subroutine create_communicator(new_comm, num_procs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ocean communicator.
!
!  this routine should be called from init_domain1 when the
!  domain configuration (e.g. nprocs_btrop) has been determined
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

!EOP
!BOC
!-----------------------------------------------------------------------

   call POP_CommCreateCommunicator(new_comm, num_procs)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
