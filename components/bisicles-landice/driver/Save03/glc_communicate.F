!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module glc_communicate

! !MODULE: glc_communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.
!
!  WJS (11-19-12): some information here is redundant with information in glimmer-cism's
!  parallel module - e.g., my_task (redundant with this_rank) and the get_num_procs
!  routine. However, I am keeping this redundant information here, looking to the future:
!  When we have multiple instances of cism, and/or GIC, all within a single GLC: the
!  information here will tell us about the MPI information relative to the whole GLC
!  communicator, whereas the information in CISM's parallel module will tell us about the
!  MPI information relative to CISM's MPI communicator, which could theoretically be a
!  sub-communicator of MPI_COMM_GLC.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_communicate.F90 66 2007-05-02 16:52:51Z dbailey $
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
!
! !USES:

   use glc_kinds_mod
   use shr_sys_mod, only : shr_sys_abort

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
      MPI_COMM_GLC,             &! MPI communicator for glc comms
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

   integer (int_kind), parameter, public :: &
      mpitag_bndy_2d        = 1,    &! MPI tags for various
      mpitag_bndy_3d        = 2,    &! communication patterns
      mpitag_gs             = 1000   ! 

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_communicate
! !INTERFACE:

 subroutine init_communicate(mpicom)

! !DESCRIPTION:
!  This routine sets up MPI environment and defines glc communicator.
!
! !REVISION HISTORY:
!  same as module
!
! !USES:
   use parallel, only : parallel_set_info

! !INPUT PARAMETERS:
 
   integer (int_kind), intent(in) :: mpicom   ! MPI communicator

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   include 'mpif.h'   ! MPI Fortran include file

   integer (int_kind) :: ierr  ! MPI error flag

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ocean communications
!
!-----------------------------------------------------------------------

   MPI_COMM_GLC = mpicom
   master_task = 0
   call MPI_COMM_RANK  (MPI_COMM_GLC, my_task, ierr)

   call parallel_set_info(MPI_COMM_GLC, master_task)

!-----------------------------------------------------------------------
!
!  On some 64-bit machines where real_kind and dbl_kind are
!  identical, the MPI implementation uses MPI_REAL for both.
!  In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------

   MPI_DBL = MPI_DOUBLE_PRECISION

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
!  MPI_COMM_GLC
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

   call MPI_COMM_SIZE(MPI_COMM_GLC, get_num_procs, ierr)

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

  return 
 
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
 
!   call MPI_BARRIER(MPI_COMM_GLC,ierr)
!   ierr = 13
!   call MPI_ABORT(0,ierr)
   call shr_sys_abort('glc_communicate.F90: abort_message_environment')
 
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
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     MPI_GROUP_GLC,         &! group of processors assigned to glc
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

   call MPI_COMM_GROUP (MPI_COMM_GLC, MPI_GROUP_GLC, ierr)

   range(1) = 0
   range(2) = num_procs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in MPI_COMM_GLC
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_GLC, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (MPI_COMM_GLC, MPI_GROUP_NEW,  &
                         new_comm, ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module glc_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
