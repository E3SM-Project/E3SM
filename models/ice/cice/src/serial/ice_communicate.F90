!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module ice_communicate

! !MODULE: ice_communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.  This instance of the module
!  is for serial execution so not much is done.
!
! !REVISION HISTORY:
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
!
! !USES:

   use ice_kinds_mod

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
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

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
!  This routine sets up MPI environment and defines ice communicator.
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

#ifdef coupled
   include 'mpif.h'   ! MPI Fortran include file

   integer (int_kind) :: ierr  ! MPI error flag
#endif

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ice communications
!
!-----------------------------------------------------------------------

#ifdef coupled
   call MPI_INIT(ierr)
   call MPI_COMM_RANK  (MPI_COMM_ICE, my_task, ierr)
#else
   my_task = 0
#endif

   master_task = 0

#ifdef coupled
!-----------------------------------------------------------------------
!
!  On some 64-bit machines where real_kind and dbl_kind are
!  identical, the MPI implementation uses MPI_REAL for both.
!  In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------

   MPI_DBL = MPI_DOUBLE_PRECISION

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: get_num_procs
! !INTERFACE:

 function get_num_procs()

! !DESCRIPTION:
!  This function returns the number of processors assigned to
!  the ice model.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  serial execution, must be only 1
!
!-----------------------------------------------------------------------

   get_num_procs = 1

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

#ifdef coupled
! !INCLUDES:

   include 'mpif.h'

#endif
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

!EOP
!BOC
#ifdef coupled
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

#else
   new_comm = MPI_COMM_ICE
#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module ice_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
