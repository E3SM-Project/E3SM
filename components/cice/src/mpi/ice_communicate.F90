!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module ice_communicate

! !MODULE: ice_communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_communicate.F90 127 2008-04-28 21:59:36Z eclare $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
!
! !USES:

   use ice_kinds_mod

#if defined key_oasis3
   use cpl_oasis3
#endif

#if defined key_oasis4
   use cpl_oasis4
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_communicate,          &
              get_num_procs,             &
              create_communicator

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      MPI_COMM_ICE,             &! MPI communicator for ice comms
      mpiR16,                   &! MPI type for r16_kind
      mpiR8,                    &! MPI type for dbl_kind
      mpiR4,                    &! MPI type for real_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

   integer (int_kind), parameter, public :: &
      mpitagHalo            = 1,    &! MPI tags for various
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
!  This routine sets up MPI environment and defines ice
!  communicator.
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

   include 'mpif.h'   ! MPI Fortran include file

   integer (int_kind) :: ierr  ! MPI error flag

   integer (int_kind) :: ice_comm

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ice communications
!
!-----------------------------------------------------------------------

#if (defined key_oasis3 || defined key_oasis4)
    ice_comm = localComm       ! communicator from NEMO/OASISn 
#else
    ice_comm = MPI_COMM_WORLD  ! Global communicator 
#endif 

#if (defined popcice || defined CICE_IN_NEMO)
   ! MPI_INIT is called elsewhere in coupled configuration
#else
   call MPI_INIT(ierr)
#endif

   call MPI_BARRIER (ice_comm, ierr)
   call MPI_COMM_DUP(ice_comm, MPI_COMM_ICE, ierr)

   master_task = 0
   call MPI_COMM_RANK  (MPI_COMM_ICE, my_task, ierr)

   mpiR16 = MPI_REAL16
   mpiR8  = MPI_REAL8
   mpiR4  = MPI_REAL4

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
!  MPI_COMM_ICE
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr

!-----------------------------------------------------------------------

   call MPI_COMM_SIZE(MPI_COMM_ICE, get_num_procs, ierr)

!-----------------------------------------------------------------------
!EOC

 end function get_num_procs

!***********************************************************************
!BOP
! !IROUTINE: create_communicator
! !INTERFACE:

 subroutine create_communicator(new_comm, num_procs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ice communicator.
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
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     MPI_GROUP_ICE,         &! group of processors assigned to ice
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (int_kind) :: &
     ierr                    ! error flag for MPI comms

   integer (int_kind), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP (MPI_COMM_ICE, MPI_GROUP_ICE, ierr)

   range(1) = 0
   range(2) = num_procs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in MPI_COMM_ICE
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_ICE, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (MPI_COMM_ICE, MPI_GROUP_NEW,  &
                         new_comm, ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module ice_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
