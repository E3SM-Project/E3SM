! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_dmpar
!
!> \brief MPAS Communication Routines
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module contains all communication routines. All MPI calls should be made in this module.
!
!-----------------------------------------------------------------------
module mpas_dmpar

#define COMMA ,
#define DMPAR_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define DMPAR_WARNING_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define DMPAR_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

   use mpas_derived_types
   use mpas_sort
   use mpas_hash
   use mpas_io_units
   use mpas_threading
   use mpas_pool_routines

#ifdef _MPI
#ifndef NOMPIMOD
   use mpi
#endif
#endif

   implicit none
   private

#ifdef _MPI
#ifdef NOMPIMOD
include 'mpif.h'
#endif
   integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER
   integer, parameter :: MPI_2INTEGERKIND = MPI_2INTEGER

#ifdef SINGLE_PRECISION
   integer, parameter :: MPI_REALKIND = MPI_REAL
   integer, parameter :: MPI_2REALKIND = MPI_2REAL
#else
   integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
   integer, parameter :: MPI_2REALKIND = MPI_2DOUBLE_PRECISION
#endif
#endif

   integer, parameter, public :: IO_NODE = 0
   integer, parameter, public :: BUFSIZE = 6000

   public :: mpas_dmpar_init
   public :: mpas_dmpar_finalize
   public :: mpas_dmpar_abort
   public :: mpas_dmpar_bcast_int
   public :: mpas_dmpar_bcast_ints
   public :: mpas_dmpar_bcast_real
   public :: mpas_dmpar_bcast_reals
   public :: mpas_dmpar_bcast_double
   public :: mpas_dmpar_bcast_doubles
   public :: mpas_dmpar_bcast_single
   public :: mpas_dmpar_bcast_singles
   public :: mpas_dmpar_bcast_logical
   public :: mpas_dmpar_bcast_char
   public :: mpas_dmpar_bcast_chars
   public :: mpas_dmpar_sum_int
   public :: mpas_dmpar_sum_real
   public :: mpas_dmpar_min_int
   public :: mpas_dmpar_min_real
   public :: mpas_dmpar_max_int
   public :: mpas_dmpar_max_real
   public :: mpas_dmpar_minloc_int
   public :: mpas_dmpar_minloc_real
   public :: mpas_dmpar_maxloc_int
   public :: mpas_dmpar_maxloc_real
   public :: mpas_dmpar_sum_int_array
   public :: mpas_dmpar_min_int_array
   public :: mpas_dmpar_max_int_array
   public :: mpas_dmpar_sum_real_array
   public :: mpas_dmpar_min_real_array
   public :: mpas_dmpar_max_real_array
   public :: mpas_dmpar_scatter_ints
   public :: mpas_dmpar_get_index_range
   public :: mpas_dmpar_compute_index_range
   public :: mpas_dmpar_get_exch_list
   public :: mpas_dmpar_build_comm_lists
   public :: mpas_dmpar_init_multihalo_exchange_list
   public :: mpas_dmpar_destroy_mulithalo_exchange_list
   public :: mpas_dmpar_destroy_communication_list
   public :: mpas_dmpar_destroy_exchange_list
   public :: mpas_dmpar_get_time
   public :: mpas_dmpar_exch_group_create
   public :: mpas_dmpar_exch_group_destroy
   public :: mpas_dmpar_exch_group_add_field
   public :: mpas_dmpar_exch_group_remove_field
   public :: mpas_dmpar_exch_group_begin_halo_exch
   public :: mpas_dmpar_exch_group_local_halo_exch
   public :: mpas_dmpar_exch_group_end_halo_exch
   public :: mpas_dmpar_exch_group_full_halo_exch
   public :: mpas_dmpar_field_halo_exch


   interface mpas_dmpar_alltoall_field
      module procedure mpas_dmpar_alltoall_field1d_integer
      module procedure mpas_dmpar_alltoall_field2d_integer
      module procedure mpas_dmpar_alltoall_field1d_real
      module procedure mpas_dmpar_alltoall_field2d_real
      module procedure mpas_dmpar_alltoall_field3d_real
      module procedure mpas_dmpar_alltoall_field4d_real
      module procedure mpas_dmpar_alltoall_field5d_real
   end interface

   public :: mpas_dmpar_alltoall_field

   private :: mpas_dmpar_alltoall_field1d_integer
   private :: mpas_dmpar_alltoall_field2d_integer
   private :: mpas_dmpar_alltoall_field1d_real
   private :: mpas_dmpar_alltoall_field2d_real
   private :: mpas_dmpar_alltoall_field3d_real
   private :: mpas_dmpar_alltoall_field4d_real
   private :: mpas_dmpar_alltoall_field5d_real


   interface mpas_dmpar_exch_halo_field
      module procedure mpas_dmpar_exch_halo_field1d_integer
      module procedure mpas_dmpar_exch_halo_field2d_integer
      module procedure mpas_dmpar_exch_halo_field3d_integer
      module procedure mpas_dmpar_exch_halo_field1d_real
      module procedure mpas_dmpar_exch_halo_field2d_real
      module procedure mpas_dmpar_exch_halo_field3d_real
      module procedure mpas_dmpar_exch_halo_field4d_real
      module procedure mpas_dmpar_exch_halo_field5d_real
   end interface

   public :: mpas_dmpar_exch_halo_field

   private :: mpas_dmpar_exch_halo_field1d_integer
   private :: mpas_dmpar_exch_halo_field2d_integer
   private :: mpas_dmpar_exch_halo_field3d_integer
   private :: mpas_dmpar_exch_halo_field1d_real
   private :: mpas_dmpar_exch_halo_field2d_real
   private :: mpas_dmpar_exch_halo_field3d_real
   private :: mpas_dmpar_exch_halo_field4d_real
   private :: mpas_dmpar_exch_halo_field5d_real

   interface mpas_dmpar_copy_field
      module procedure mpas_dmpar_copy_field1d_integer
      module procedure mpas_dmpar_copy_field2d_integer
      module procedure mpas_dmpar_copy_field3d_integer
      module procedure mpas_dmpar_copy_field1d_real
      module procedure mpas_dmpar_copy_field2d_real
      module procedure mpas_dmpar_copy_field3d_real
      module procedure mpas_dmpar_copy_field4d_real
      module procedure mpas_dmpar_copy_field5d_real
   end interface

   public :: mpas_dmpar_copy_field

   private :: mpas_dmpar_copy_field1d_integer
   private :: mpas_dmpar_copy_field2d_integer
   private :: mpas_dmpar_copy_field3d_integer
   private :: mpas_dmpar_copy_field1d_real
   private :: mpas_dmpar_copy_field2d_real
   private :: mpas_dmpar_copy_field3d_real
   private :: mpas_dmpar_copy_field4d_real
   private :: mpas_dmpar_copy_field5d_real

   interface mpas_dmpar_exch_group_pack_buffer_field
      module procedure mpas_dmpar_exch_group_pack_buffer_field1d_integer
      module procedure mpas_dmpar_exch_group_pack_buffer_field2d_integer
      module procedure mpas_dmpar_exch_group_pack_buffer_field3d_integer
      module procedure mpas_dmpar_exch_group_pack_buffer_field1d_real
      module procedure mpas_dmpar_exch_group_pack_buffer_field2d_real
      module procedure mpas_dmpar_exch_group_pack_buffer_field3d_real
      module procedure mpas_dmpar_exch_group_pack_buffer_field4d_real
      module procedure mpas_dmpar_exch_group_pack_buffer_field5d_real
   end interface

   interface mpas_dmpar_exch_group_local_exch_field
      module procedure mpas_dmpar_exch_group_local_exch_field1d_integer
      module procedure mpas_dmpar_exch_group_local_exch_field2d_integer
      module procedure mpas_dmpar_exch_group_local_exch_field3d_integer
      module procedure mpas_dmpar_exch_group_local_exch_field1d_real
      module procedure mpas_dmpar_exch_group_local_exch_field2d_real
      module procedure mpas_dmpar_exch_group_local_exch_field3d_real
      module procedure mpas_dmpar_exch_group_local_exch_field4d_real
      module procedure mpas_dmpar_exch_group_local_exch_field5d_real
   end interface mpas_dmpar_exch_group_local_exch_field

   interface mpas_dmpar_exch_group_unpack_buffer_field
      module procedure mpas_dmpar_exch_group_unpack_buffer_field1d_integer
      module procedure mpas_dmpar_exch_group_unpack_buffer_field2d_integer
      module procedure mpas_dmpar_exch_group_unpack_buffer_field3d_integer
      module procedure mpas_dmpar_exch_group_unpack_buffer_field1d_real
      module procedure mpas_dmpar_exch_group_unpack_buffer_field2d_real
      module procedure mpas_dmpar_exch_group_unpack_buffer_field3d_real
      module procedure mpas_dmpar_exch_group_unpack_buffer_field4d_real
      module procedure mpas_dmpar_exch_group_unpack_buffer_field5d_real
   end interface



   contains

!-----------------------------------------------------------------------
!  routine mpas_dmpar_init
!
!> \brief MPAS dmpar initialization routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine initializes dmpar. It calls MPI_Init (if required), and setups up the communicators.
!>  It also setups of the domain information structure.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_init(dminfo, mpi_comm)!{{{

      implicit none

      type (dm_info), intent(inout) :: dminfo !< Input/Output: Domain information
      integer, intent(in), optional :: mpi_comm !< Input - Optional: externally-supplied MPI communicator

#ifdef _MPI
      integer :: mpi_rank, mpi_size
      integer :: mpi_ierr
#ifdef MPAS_OPENMP
      integer :: desiredThreadLevel, threadLevel
#endif

      if (present(mpi_comm)) then
         dminfo % comm = mpi_comm
         dminfo % using_external_comm = .true.
      else
#ifdef MPAS_OPENMP
         desiredThreadLevel = MPI_THREAD_MULTIPLE
         call MPI_Init_thread(desiredThreadLevel, threadLevel, mpi_ierr)
         if ( threadLevel /= desiredThreadLevel ) then
            write(stderrUnit, *) 'Warning: MPI implementation gave thread level of ', threadLevel, ' when ', desiredThreadLevel, ' was requested.'
         end if
#else
         call MPI_Init(mpi_ierr)
#endif
         dminfo % comm = MPI_COMM_WORLD
         dminfo % using_external_comm = .false.
      end if

      ! Find out our rank and the total number of processors
      call MPI_Comm_rank(dminfo % comm, mpi_rank, mpi_ierr)
      call MPI_Comm_size(dminfo % comm, mpi_size, mpi_ierr)

      dminfo % nprocs = mpi_size
      dminfo % my_proc_id = mpi_rank

      write(stderrUnit,'(a,i5,a,i5,a)') 'task ', mpi_rank, ' of ', mpi_size, &
        ' is running'

      call open_streams(dminfo % my_proc_id)

      dminfo % info = MPI_INFO_NULL
#else
      dminfo % comm = 0
      dminfo % my_proc_id = IO_NODE
      dminfo % nprocs = 1
      dminfo % using_external_comm = .false.
#endif

   end subroutine mpas_dmpar_init!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_finalize
!
!> \brief MPAS dmpar finalization routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine finalizes dmpar. It calls MPI_Finalize (if required).
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_finalize(dminfo)!{{{

      implicit none

      type (dm_info), intent(inout) :: dminfo !< Input/Output: Domain information.

#ifdef _MPI
      integer :: mpi_ierr

      if (.not. dminfo % using_external_comm) then
         call MPI_Finalize(mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_finalize!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_abort
!
!> \brief MPAS dmpar abort routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine aborts MPI. A call to it kills the model through the use of MPI_Abort.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_abort(dminfo)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      call MPI_Abort(dminfo % comm, mpi_errcode, mpi_ierr)
#endif

      stop

   end subroutine mpas_dmpar_abort!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_int
!
!> \brief MPAS dmpar broadcast integer routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts an integer to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_int(dminfo, i, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(inout) :: i !< Input/Output: Integer to broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(i, 1, MPI_INTEGERKIND, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_ints
!
!> \brief MPAS dmpar broadcast integers routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts an array of integers to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_ints(dminfo, n, iarray, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      integer, dimension(n), intent(inout) :: iarray !< Input/Output: Array of integers
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(iarray, n, MPI_INTEGERKIND, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_ints!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_real
!
!> \brief MPAS dmpar broadcast real routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts a real to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_real(dminfo, r, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real (kind=RKIND), intent(inout) :: r !< Input/Output: Real to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(r, 1, MPI_REALKIND, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_reals
!
!> \brief MPAS dmpar broadcast reals routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts an array of reals to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_reals(dminfo, n, rarray, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      real (kind=RKIND), dimension(n), intent(inout) :: rarray !< Input/Output: Array of reals to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(rarray, n, MPI_REALKIND, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_reals!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_double
!
!> \brief MPAS dmpar broadcast double routine.
!> \author Michael Duda
!> \date   11/04/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts a double to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_double(dminfo, r, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      double precision, intent(inout) :: r !< Input/Output: Double to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(r, 1, MPI_DOUBLE_PRECISION, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_double!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_doubles
!
!> \brief MPAS dmpar broadcast doubles routine.
!> \author Michael Duda
!> \date   11/04/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts an array of doubles to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_doubles(dminfo, n, rarray, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      double precision, dimension(n), intent(inout) :: rarray !< Input/Output: Array of doubles to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(rarray, n, MPI_DOUBLE_PRECISION, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_doubles!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_double
!
!> \brief MPAS dmpar broadcast double routine.
!> \author Michael Duda
!> \date   11/04/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts a double to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_single(dminfo, r, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=R4KIND), intent(inout) :: r !< Input/Output: Double to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(r, 1, MPI_REAL4, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_single!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_doubles
!
!> \brief MPAS dmpar broadcast doubles routine.
!> \author Michael Duda
!> \date   11/04/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts an array of doubles to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_singles(dminfo, n, rarray, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      real(kind=R4KIND), dimension(n), intent(inout) :: rarray !< Input/Output: Array of doubles to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(rarray, n, MPI_REAL4, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_singles!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_logical
!
!> \brief MPAS dmpar broadcast logical routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts a logical to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_logical(dminfo, l, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      logical, intent(inout) :: l !< Input/Output: Logical to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: itemp
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         if (dminfo % my_proc_id == IO_NODE) then
            if (l) then
               itemp = 1
            else
               itemp = 0
            end if
         end if

         call MPI_Bcast(itemp, 1, MPI_INTEGERKIND, source, dminfo % comm, mpi_ierr)

         if (itemp == 1) then
            l = .true.
         else
            l = .false.
         end if
      end if
#endif

   end subroutine mpas_dmpar_bcast_logical!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_char
!
!> \brief MPAS dmpar broadcast character routine.
!> \author Michael Duda
!> \date   03/26/13; modified by William Lipscomb 01/21/15
!> \details
!>  This routine broadcasts a character to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_char(dminfo, c, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      character (len=*), intent(inout) :: c !< Input/Output: Character to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(c, len(c), MPI_CHARACTER, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_char!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_chars
!
!> \brief MPAS dmpar broadcast character array routine.
!> \author Doug Jacobsen
!> \date   01/22/2016
!> \details
!>  This routine broadcasts an array of characters to all processors in the communicator.
!>  An optional argument specifies the source node; else broadcast from IO_NODE.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_chars(dminfo, n, carray, proc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Number of character strings to broadcast
      character (len=*), dimension(:), intent(inout) :: carray !< Input/Output: Character to be broadcast
      integer, intent(in), optional :: proc  !< optional argument indicating which processor to broadcast from

#ifdef _MPI
      integer :: mpi_ierr, source
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(proc)) then
            source = proc
         else
            source = IO_NODE
         endif

         call MPI_Bcast(carray, n * len(carray(1)), MPI_CHARACTER, source, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_bcast_chars!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_sum_int
!
!> \brief MPAS dmpar sum integers routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine sums (Allreduce) integer values across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_sum_int(dminfo, i, isum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: i !< Input: Integer value input
      integer, intent(out) :: isum !< Output: Integer sum for output

      integer :: mpi_ierr
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(i, isum, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
         isum = i
#endif
      end if

   end subroutine mpas_dmpar_sum_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_sum_real
!
!> \brief MPAS dmpar sum real routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine sums (Allreduce) real values across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_sum_real(dminfo, r, rsum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=RKIND), intent(in) :: r !< Input: Real values to be summed
      real(kind=RKIND), intent(out) :: rsum  !< Output: Sum of reals for output

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(r, rsum, 1, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
         rsum = r
#endif
      end if

   end subroutine mpas_dmpar_sum_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_min_int
!
!> \brief MPAS dmpar minimum integer routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine returns the minimum integer value across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_min_int(dminfo, i, imin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: i !< Input: Integer value
      integer, intent(out) :: imin !< Output: Minimum integer value

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(i, imin, 1, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
         imin = i
#endif
      end if

   end subroutine mpas_dmpar_min_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_min_real
!
!> \brief MPAS dmpar minimum real routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine returns the minimum real value across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_min_real(dminfo, r, rmin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=RKIND), intent(in) :: r !< Input: Real value
      real(kind=RKIND), intent(out) :: rmin !< Output: Minimum of real value

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(r, rmin, 1, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
         rmin = r
#endif
      end if

   end subroutine mpas_dmpar_min_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_max_int
!
!> \brief MPAS dmpar maximum integer routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine returns the maximum integer value across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_max_int(dminfo, i, imax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: i !< Input: Integer value
      integer, intent(out) :: imax !< Output: Maximum of integer values

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(i, imax, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
         imax = i
#endif
      end if

   end subroutine mpas_dmpar_max_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_max_real
!
!> \brief MPAS dmpar maximum real routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine returns the maximum real value across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_max_real(dminfo, r, rmax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=RKIND), intent(in) :: r !< Input: Real value
      real(kind=RKIND), intent(out) :: rmax !< Output: Maximum of real values

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(r, rmax, 1, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
         rmax = r
#endif
      end if

   end subroutine mpas_dmpar_max_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_minloc_int
!
!> \brief MPAS dmpar minloc integer routine.
!> \author William Lipscomb
!> \date   01/21/15
!> \details
!>  This routine returns the minimum integer value across all processors in a communicator,
!>  along with the processor on which this value resides.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_minloc_int(dminfo, i, imin, procout)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: i !< Input: Integer value
      integer, intent(out) :: imin !< Output: Minimum of integer values
      integer, intent(out) :: procout  !< Output: Processor on which imin resides
      integer :: mpi_ierr, threadNum
      integer, dimension(2,1) :: recvbuf, sendbuf

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         sendbuf(1,1) = i
         sendbuf(2,1) = dminfo % my_proc_id  ! This is the processor number associated with the value i
         call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_2INTEGERKIND, MPI_MINLOC, dminfo % comm, mpi_ierr)
         imin = recvbuf(1,1)
         procout = recvbuf(2,1)
#else
         imin = i
         procout = IO_NODE
#endif
      end if

   end subroutine mpas_dmpar_minloc_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_minloc_real
!
!> \brief MPAS dmpar minloc real routine.
!> \author William Lipscomb
!> \date   01/21/15
!> \details
!>  This routine returns the minimum real value across all processors in a communicator,
!>  along with the processor on which this value resides.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_minloc_real(dminfo, r, rmin, procout)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=RKIND), intent(in) :: r !< Input: Real value
      real(kind=RKIND), intent(out) :: rmin !< Output: Minimum of real values
      integer, intent(out) :: procout  !< Output: Processor on which rin resides
      integer :: mpi_ierr, threadNum
      real(kind=RKIND), dimension(2,1) :: recvbuf, sendbuf

      if ( threadNum == 0 ) then
#ifdef _MPI
         sendbuf(1,1) = r
         sendbuf(2,1) = dminfo % my_proc_id  ! This is the processor number associated with the value x (coerced to a real)
         call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_2REALKIND, MPI_MINLOC, dminfo % comm, mpi_ierr)
         rmin = recvbuf(1,1)
         procout = recvbuf(2,1)   ! coerced back to integer
#else
         rmin = r
         procout = IO_NODE
#endif
      end if

   end subroutine mpas_dmpar_minloc_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_maxloc_int
!
!> \brief MPAS dmpar maxloc integer routine.
!> \author William Lipscomb
!> \date   01/21/15
!> \details
!>  This routine returns the maximum integer value across all processors in a communicator,
!>  along with the processor on which this value resides.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_maxloc_int(dminfo, i, imax, procout)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: i !< Input: Integer value
      integer, intent(out) :: imax !< Output: Maximum of integer values
      integer, intent(out) :: procout  !< Output: Processor on which imax resides
      integer :: mpi_ierr, threadNum
      integer, dimension(2,1) :: recvbuf, sendbuf

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         sendbuf(1,1) = i
         sendbuf(2,1) = dminfo % my_proc_id  ! This is the processor number associated with the value i
         call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_2INTEGERKIND, MPI_MAXLOC, dminfo % comm, mpi_ierr)
         imax = recvbuf(1,1)
         procout = recvbuf(2,1)
#else
         imax = i
         procout = IO_NODE
#endif
      end if

   end subroutine mpas_dmpar_maxloc_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_maxloc_real
!
!> \brief MPAS dmpar maxloc real routine.
!> \author William Lipscomb
!> \date   01/21/15
!> \details
!>  This routine returns the maximum real value across all processors in a communicator,
!>  along with the processor on which this value resides.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_maxloc_real(dminfo, r, rmax, procout)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real(kind=RKIND), intent(in) :: r !< Input: Real value
      real(kind=RKIND), intent(out) :: rmax !< Output: Maximum of real values
      integer, intent(out) :: procout  !< Output: Processor on which rmax resides
      integer :: mpi_ierr, threadNum
      real(kind=RKIND), dimension(2,1) :: recvbuf, sendbuf

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         sendbuf(1,1) = r
         sendbuf(2,1) = dminfo % my_proc_id  ! This is the processor number associated with the value x (coerced to a real)
         call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_2REALKIND, MPI_MAXLOC, dminfo % comm, mpi_ierr)
         rmax = recvbuf(1,1)
         procout = recvbuf(2,1)   ! coerced back to integer
#else
         rmax = r
         procout = IO_NODE
#endif
      end if

   end subroutine mpas_dmpar_maxloc_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_sum_int_array
!
!> \brief MPAS dmpar integer array sum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes the sum of a set of integer arrays across all processors in a communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_sum_int_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      integer, dimension(nElements), intent(in) :: inArray !< Input: Processor specific array to sum
      integer, dimension(nElements), intent(out) :: outArray !< Output: Sum of arrays

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_sum_int_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_min_int_array
!
!> \brief MPAS dmpar integer array minimum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes an array of minimum values for each index across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_min_int_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Array size
      integer, dimension(nElements), intent(in) :: inArray !< Input: Input array of integers
      integer, dimension(nElements), intent(out) :: outArray !< Output: Array of minimum integers

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_min_int_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_min_int_array
!
!> \brief MPAS dmpar integer array maximum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes an array of maximum values for each index across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_max_int_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      integer, dimension(nElements), intent(in) :: inArray !< Input: Array of integers
      integer, dimension(nElements), intent(out) :: outArray !< Output: Array of maximum integers

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_max_int_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_sum_real_array
!
!> \brief MPAS dmpar real array sum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes the sum array of real values  across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_sum_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray !< Input: Array of reals
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray !< Output: Array of real sums

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_sum_real_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_min_real_array
!
!> \brief MPAS dmpar real array minimum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes the minimum array of real values  across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_min_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray !< Input: Array of reals
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray !< Input: Array of minimum reals

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_min_real_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_max_real_array
!
!> \brief MPAS dmpar real array maximum routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes the maximum array of real values  across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_max_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray !< Input: Array of reals
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray !< Output: Array of maximum reals

      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
#ifdef _MPI
         call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
         outArray = inArray
#endif
      end if

   end subroutine mpas_dmpar_max_real_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_maxloc_real_array
!
!> \brief MPAS dmpar location of maximum real value
!> \author Michael Duda
!> \date   03/13/14
!> \details
!>  Returns the maximum real among arrays of reals provided by all tasks, as
!   well as the corresponding location of that maximum. If the maximum value
!   occurs in more than one location, the returned location may be any of these.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_maxloc_real_array(dminfo, nElements, inArray, inLoc, maxVal, maxLoc)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nElements !< Input: Length of arrays
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray   !< Input: Array of reals
      integer, dimension(nElements), intent(in) :: inLoc              !< Input: Array of locations
      real(kind=RKIND), intent(out) :: maxVal                         !< Output: Maximum of all reals
      integer, intent(out) :: maxLoc                                  !< Output: Maximum of all reals

      integer :: mpi_ierr
      real(kind=RKIND), dimension(2) :: in, out
      integer :: i

      !
      ! Find local max and maxloc
      ! NB: We will have a problem on our hands if the maxloc value exceeds the
      ! range of exactly representable integral values for the realkind being used
      !
      in(1) = inArray(1)
      in(2) = 1.0
      do i=2,nElements
         if (inArray(i) > in(1)) then
            in(1) = inArray(i)
            in(2) = i
         end if
      end do

#ifdef _MPI
#ifdef SINGLE_PRECISION
      call MPI_Allreduce(in, out, 1, MPI_2REAL, MPI_MAXLOC, dminfo % comm, mpi_ierr)
#else
      call MPI_Allreduce(in, out, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, dminfo % comm, mpi_ierr)
#endif
#endif
      maxVal = out(1)
      maxLoc = out(2)

   end subroutine mpas_dmpar_maxloc_real_array!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_scatter_ints
!
!> \brief MPAS dmpar scatter integers routine
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine computes the maximum array of real values  across all processors in a communicator, from some input arrays.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_scatter_ints(dminfo, nprocs, noutlist, displs, counts, inlist, outlist)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: nprocs !< Input: Number of processors
      integer, intent(in) :: noutlist !< Input: Number integers to receive
      integer, dimension(nprocs), intent(in) :: displs !< Input: Displacement in sending array
      integer, dimension(nprocs), intent(in) :: counts !< Input: Number of integers to distribute
      integer, dimension(:), pointer :: inlist !< Input: List of integers to send
      integer, dimension(noutlist), intent(inout) :: outlist !< Output: List of received integers

#ifdef _MPI
      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         call MPI_Scatterv(inlist, counts, displs, MPI_INTEGERKIND, outlist, noutlist, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_scatter_ints!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_get_index_range
!
!> \brief MPAS dmpar processor specific range of indices
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine divides a global range of indices among all processors, and returns the range of indices a specific processors is responsible for.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_get_index_range(dminfo, &!{{{
                                    global_start, global_end, &
                                    local_start, local_end)

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: global_start !< Input: Starting index in global range
      integer, intent(in) :: global_end !< Input: Ending index in global range
      integer, intent(out) :: local_start !< Output: Starting index in local range
      integer, intent(out) :: local_end !< Output: Ending index in local range

      local_start = nint(real(dminfo % my_proc_id,R8KIND) * real(global_end - global_start + 1,R8KIND) &
                       / real(dminfo % nprocs,R8KIND)) + 1
      local_end   = nint(real(dminfo % my_proc_id + 1,R8KIND) * real(global_end - global_start + 1,R8KIND) &
                       / real(dminfo % nprocs,R8KIND))

   end subroutine mpas_dmpar_get_index_range!}}}

   subroutine mpas_dmpar_compute_index_range(dminfo, &!{{{
                                        local_start, local_end, &
                                        global_start, global_end)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: local_start, local_end
      integer, intent(inout) :: global_start, global_end

      integer :: n
      integer :: mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()

      n = local_end - local_start + 1

      if (dminfo % my_proc_id == 0) then
         global_start = 1
         global_end = global_start + n - 1

#ifdef _MPI
      else if (dminfo % my_proc_id == dminfo % nprocs - 1) then
         if ( threadNum == 0 ) then
            call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()
         global_end = global_start + n - 1

      else
         if ( threadNum == 0 ) then
            call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()
         global_end = global_start + n
         if ( threadNum == 0 ) then
            call MPI_Send(global_end, 1, MPI_INTEGERKIND, dminfo % my_proc_id + 1, 0, dminfo % comm, mpi_ierr)
         end if
         call mpas_threading_barrier()
         global_end = global_end - 1
#endif

      end if


   end subroutine mpas_dmpar_compute_index_range!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_get_exch_list
!
!> \brief MPAS dmpar exchange list builder
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine builds exchange lists to communicated between the lists of owned and needed fields, over a given number of halos.
!>  Exchange lists are built into the input fields.
!
!-----------------------------------------------------------------------
subroutine mpas_dmpar_get_exch_list(haloLayer, ownedListField, neededListField, offsetListField, ownedLimitField)!{{{

      implicit none

      integer, intent(in) :: haloLayer !< Input: Halo layer to build exchange list for
      type (field1dInteger), pointer :: ownedListField !< Input/Output: List of owned fields
      type (field1dInteger), pointer :: neededListField !< Input/Output: List of needed fields
      type (field0dInteger), pointer, optional :: offsetListField !< Input: Offsets for placement of received data into destination arrays
      type (field0dInteger), pointer, optional :: ownedLimitField !< Input: List of limits in owned array

      type (dm_info), pointer :: dminfo

      integer :: i, j, k, kk, iBlock
      integer :: totalSize, nMesgRecv, nMesgSend, recvNeighbor, sendNeighbor, currentProc
      integer :: totalSent
      integer, allocatable, dimension(:) :: numToSend, numToRecv
      integer, allocatable, dimension(:) :: ownedList, ownedListIndex, ownedBlock, neededList, neededListIndex, neededBlock
      integer, allocatable, dimension(:) :: offsetList, ownedLimitList
      integer, allocatable, dimension(:,:) :: ownedListSorted, ownedBlockSorted, recipientList
      integer, allocatable, dimension(:) :: ownerListIn, ownerListOut
      integer, allocatable, dimension(:) :: packingOrder
      type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
      type (field1dInteger), pointer :: fieldCursor, fieldCursor2
      type (field0dInteger), pointer :: offsetCursor, ownedLimitCursor
      integer :: nOwnedBlocks, nNeededBlocks
      integer :: nOwnedList, nNeededList
      integer :: mpi_ierr, mpi_rreq, mpi_sreq

      type (hashtable) :: neededHash
      integer :: nUniqueNeededList, threadNum
      integer, dimension(:,:), pointer :: uniqueSortedNeededList


      !
      ! *** NB: This code assumes that block % blockID values are local block IDs and are in the range [1, numBlocks]
      !         where numBlocks is the number of blocks owned by each task
      !


      ! For the ownedListField:
      !    - ownedList contains a list of the global indices owned by all blocks
      !    - ownedListIndex contains a list of the block-local indices of the global indices owned by all blocks
      !    - ownedBlock contains the local block ID associated with each index
      !
      ! Example:
      !    ownedList      := ( 21 13 15 01 05 06 33 42 44 45 )     ! Global indices from all blocks on this task
      !    ownedListIndex := (  1  2  3  4  1  2  3  4  5  6 )     ! Local  indices of global indices on each block
      !    ownedBlock     := (  1  1  1  1  2  2  2  2  2  2 )     ! Local  indices of global indices on each block
      !

      ! For the neededListField:
      !    similar to the ownedListField...

      dminfo => ownedListField % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      !
      ! Determine total number of owned blocks on this task
      !
      if ( threadNum == 0 ) then
        nOwnedBlocks = 0
        fieldCursor => ownedListField
        do while (associated(fieldCursor))
          nOwnedBlocks = nOwnedBlocks + 1
          if(associated(fieldCursor % sendList % halos(haloLayer) % exchList)) then
            call mpas_dmpar_destroy_exchange_list(fieldCursor % sendList % halos(haloLayer) % exchList)
          end if

          if(associated(fieldCursor % copyList % halos(haloLayer) % exchList)) then
            call mpas_dmpar_destroy_exchange_list(fieldCursor % copyList % halos(haloLayer) % exchList)
          end if
          fieldCursor => fieldCursor % next
        end do

        !
        ! Determine total number of needed indices on this task
        !
        nNeededList = 0
        nNeededBlocks = 0
        fieldCursor => neededListField
        do while (associated(fieldCursor))
          nNeededBlocks = nNeededBlocks + 1
          nNeededList = nNeededList + fieldCursor % dimSizes(1)
          if(associated(fieldCursor % recvList % halos(haloLayer) % exchList)) then
            call mpas_dmpar_destroy_exchange_list(fieldCursor % recvList % halos(haloLayer) % exchList)
          end if

          fieldCursor => fieldCursor % next
        end do

        !
        ! Determine unique list of needed elements.
        !
        nUniqueNeededList = 0
        call mpas_hash_init(neededHash)
        fieldCursor => neededListField
        do while (associated(fieldCursor))
          do i = 1, fieldCursor % dimSizes(1)
            if(.not. mpas_hash_search(neededHash, fieldCursor % array(i))) then
              nUniqueNeededList = nUniqueNeededList + 1
              call mpas_hash_insert(neededHash, fieldCursor % array(i))
            end if
          end do
          fieldCursor => fieldCursor % next
        end do

        kk = mpas_hash_size(neededHash)

        nUniqueNeededList = mpas_hash_size(neededHash)
        allocate(uniqueSortedNeededList(2,nUniqueNeededList))
        allocate(packingOrder(nUniqueNeededList))
        call mpas_hash_destroy(neededHash)
        call mpas_hash_init(neededHash)

        j = 0
        fieldCursor => neededListField
        do while (associated(fieldCursor))
          do i = 1, fieldCursor % dimSizes(1)
            if(.not. mpas_hash_search(neededHash, fieldCursor % array(i))) then
              j = j +1
              uniqueSortedNeededList(1, j) = fieldCursor % array(i)
              uniqueSortedNeededList(2, j) = fieldCursor % block % localBlockID
              call mpas_hash_insert(neededHash, fieldCursor % array(i))
            end if
          end do
          fieldCursor => fieldCursor % next
        end do

        kk = mpas_hash_size(neededHash)

        call mpas_hash_destroy(neededHash)
        call mpas_quicksort(nUniqueNeededList, uniqueSortedNeededList)

        !
        ! Get list of index offsets for all blocks
        !
        allocate(offsetList(nNeededBlocks))
        if (present(offsetListField)) then
          offsetCursor => offsetListField
          do while (associated(offsetCursor))
            offsetList(offsetCursor % block % localBlockID+1) = offsetCursor % scalar
            offsetCursor => offsetCursor % next
          end do
        else
          offsetList(:) = 0
        end if

        !
        ! Get list of bounds limit for owned elements
        !
        allocate(ownedLimitList(nOwnedBlocks))
        if(present(ownedLimitField)) then
          ownedLimitCursor => ownedLimitField
          do while(associated(ownedLimitCursor))
            ownedLimitList(ownedLimitCursor % block % localBlockID+1) = ownedLimitCursor % scalar
            ownedLimitCursor => ownedLimitCursor % next
          end do
        else
          fieldCursor => ownedListField
          do while(associated(fieldCursor))
            ownedLimitList(fieldCursor % block % localBlockID+1) = fieldCursor % dimSizes(1)
            fieldCursor => fieldCursor % next
          end do
        end if

        !
        ! Determine total number of owned indices on this task, and
        !   initialize output send and recv lists for ownedListField
        !
        nOwnedList = 0
        fieldCursor => ownedListField
        do while (associated(fieldCursor))
          iBlock = fieldcursor % block % localBlockID + 1
          nOwnedList = nOwnedList + ownedLimitList(iBlock)
          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI
        !
        ! Gather list of all owned indices and their associated blocks on this task
        !
        allocate(ownedList(nOwnedList))
        allocate(ownedBlock(nOwnedList))
        ownedBlock = -1
        ownedList = -1
        fieldCursor => ownedListField
        i = 1
        do while (associated(fieldCursor))
          iBlock = fieldCursor % block % localBlockID + 1
          ownedList(i:i+ownedLimitList(iBlock)-1) = fieldCursor % array(1:ownedLimitList(iBlock))
          ownedBlock(i:i+ownedLimitList(iBlock)-1) = fieldCursor % block % localBlockID
          i = i + ownedLimitList(iBlock)
          fieldCursor => fieldCursor % next
        end do

        !
        ! Gather list of all needed indices and their associated blocks on this task
        !
        allocate(neededList(nNeededList))
        allocate(neededBlock(nNeededList))
        fieldCursor => neededListField
        i = 1
        do while (associated(fieldCursor))
          neededList(i:i+fieldCursor % dimSizes(1)-1) = fieldCursor % array(:)
          neededBlock(i:i+fieldCursor % dimSizes(1)-1) = fieldCursor % block % localBlockID
          i = i + fieldCursor % dimSizes(1)
          fieldCursor => fieldCursor % next
        end do

        !
        ! Obtain sorted list of global indices owned by this task and the associated local indices and block IDs
        !
        allocate(ownedListIndex(nOwnedList))
        allocate(ownedListSorted(2,nOwnedList))
        allocate(recipientList(2,nOwnedList))
        j = 1
        k = 1
        do i=1,nOwnedList
          ownedListSorted(1,i) = ownedList(i)
          if (i > 1) then
            if(ownedBlock(i) /= ownedBlock(i-1)) k = 1
          end if
          ownedListIndex(i) = k
          ownedListSorted(2,i) = j
          j = j + 1
          k = k + 1
        end do
        call mpas_quicksort(nOwnedList, ownedListSorted)

        allocate(ownedBlockSorted(2,nOwnedList))
        do i=1,nOwnedList
          ownedBlockSorted(1,i) = ownedList(i)
          ownedBlockSorted(2,i) = ownedBlock(i)
        end do
        call mpas_quicksort(nOwnedList, ownedBlockSorted)


        allocate(neededListIndex(nNeededList))
        j = 1
        do i=1,nNeededList
          if (i > 1) then
            if(neededBlock(i) /= neededBlock(i-1)) j = 1
          end if
          neededListIndex(i) = j
          j = j + 1
        end do

        !
        ! Set totalSize to the maximum number of items in any task's needed list
        !
        call MPI_Allreduce(nUniqueNeededList, totalSize, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)

        allocate(ownerListIn(totalSize))
        allocate(ownerListOut(totalSize))

        nMesgSend = nUniqueNeededList
        nMesgRecv = nUniqueNeededList
        ownerListOut(1:nUniqueNeededList) = uniqueSortedNeededList(1,1:nUniqueNeededList)

        recvNeighbor = mod(dminfo % my_proc_id + dminfo % nprocs - 1, dminfo % nprocs)
        sendNeighbor = mod(dminfo % my_proc_id + 1, dminfo % nprocs)

        allocate(numToSend(nOwnedBlocks))
        allocate(numToRecv(nNeededBlocks))

        ! Initial send of data to neighbors.
        if(dminfo % nProcs == 1) then
          ownerListIn = ownerListOut
        else
          call MPI_Irecv(nMesgRecv, 1, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
          call MPI_Isend(nMesgSend, 1, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
          call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
          call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
          call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
        end if

        !
        ! For each processor (not including ourself), mark the indices that we will provide to
        !    that processor in ownerListOut, and build a send list for that processor if we
        !    do need to send any indices
        !
        do i=2, dminfo % nprocs
          recipientList = -1
          numToSend(:) = 0
          totalSent = 0

          currentProc = mod(dminfo % my_proc_id + dminfo % nprocs - i + 1, dminfo % nprocs)
          do j=1,nMesgRecv
            if (ownerListIn(j) > 0) then
              k = mpas_binary_search(ownedListSorted, 2, 1, nOwnedList, ownerListIn(j))
              if (k <= nOwnedList) then
                iBlock = ownedBlock(ownedListSorted(2,k)) + 1
                numToSend(iBlock) = numToSend(iBlock) + 1
                totalSent = totalSent + 1

                ! recipientList(1,:) represents the index in the srcList to place this data
                recipientList(1,ownedListSorted(2,k)) = numToSend(iBlock)
                ! recipientList(2,:) represnets the index in the buffer to place this data
                recipientList(2,ownedListSorted(2,k)) = totalSent

                ownerListOut(j) = -1 * dminfo % my_proc_id
              else
                ownerListOut(j) = ownerListIn(j)
              end if
            else
              ownerListOut(j) = ownerListIn(j)
            end if
          end do

          fieldCursor => ownedListField
          do while (associated(fieldCursor))
            iBlock = fieldCursor % block % localBlockID + 1

            if (numToSend(iBlock) > 0) then
              ! Find end of send list
              if(.not.associated(fieldCursor % sendList % halos(haloLayer) % exchList)) then
                allocate(fieldCursor % sendList % halos(haloLayer) % exchList)
                exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
                nullify(exchListPtr % next)
              else
                exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
                exchListPtr2 => fieldCursor % sendList % halos(haloLayer) % exchList % next
                do while(associated(exchListPtr2))
                  exchListPtr => exchListPtr % next
                  exchListPtr2 => exchListPtr % next
                end do

                allocate(exchListPtr % next)
                exchListPtr => exchListPtr % next
                nullify(exchListPtr % next)
              end if

              exchListPtr % endPointID = currentProc
              exchListPtr % nlist = numToSend(iBlock)
              allocate(exchListPtr % srcList(numToSend(iBlock)))
              allocate(exchListPtr % destList(numToSend(iBlock)))
              exchListPtr % srcList = -1
              exchListPtr % destList = -1

              kk = 1
              do j=1,nOwnedList
                if (recipientList(1,j) /= -1) then
                  if(ownedBlock(j) == fieldCursor % block % localBlockID) then
                    exchListPtr % srcList(recipientList(1,j)) = ownedListIndex(j)
                    exchListPtr % destList(recipientList(1,j)) = recipientList(2,j)
                    kk = kk + 1
                  end if
                end if
              end do
            end if

            fieldCursor => fieldCursor % next
          end do

          nMesgSend = nMesgRecv
          call MPI_Irecv(nMesgRecv, 1, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
          call MPI_Isend(nMesgSend, 1, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
          call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
          call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
          call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
          call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
        end do

        !
        ! With our needed list returned to us, build receive lists based on which indices were
        !    marked by other tasks
        !
        do i=0, dminfo % nprocs - 1
          if(i == dminfo % my_proc_id) cycle

          numToRecv(:) = 0
          packingOrder = 0

          k = 0
          do j=1,nUniqueNeededList
            if (ownerListIn(j) == -i) then
              k = k + 1
              packingOrder(j) = k
            end if
          end do

          fieldCursor => neededListField
          do while (associated(fieldCursor))
            do j = 1, fieldCursor % dimSizes(1)
              k = mpas_binary_search(uniqueSortedNeededList, 2, 1, nUniqueNeededList, fieldCursor % array(j))
              if(k <= nUniqueNeededList) then
                if(ownerListIn(k) == -i) then
                  iBlock = fieldCursor % block % localBlockID + 1
                  numToRecv(iBlock) = numToRecv(iBlock) + 1
                end if
              end if
            end do
            fieldCursor => fieldCursor % next
          end do

          fieldCursor => neededListField
          do while (associated(fieldCursor))
            iBlock = fieldCursor % block % localBlockID + 1

            if (numToRecv(iBlock) > 0) then
              if(.not.associated(fieldCursor % recvList % halos(haloLayer) % exchList)) then
                allocate(fieldCursor % recvList % halos(haloLayer) % exchList)
                exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
                nullify(exchListPtr % next)
              else
                ! Find end of recv list
                exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
                exchListPtr2 => fieldCursor % recvList % halos(haloLayer) % exchList % next
                do while(associated(exchListPtr2))
                  exchListPtr => exchListPtr % next
                  exchListPtr2 => exchListPtr % next
                end do

                allocate(exchListPtr % next)
                exchListPtr => exchListPtr % next
                nullify(exchListPtr % next)
              end if

              exchListPtr % endPointID = i
              exchListPtr % nlist = numToRecv(iBlock)
              allocate(exchListPtr % srcList(exchListPtr % nList))
              allocate(exchListPtr % destList(exchListPtr % nList))
              exchListPtr % srcList = -1
              exchListPtr % destList = -1

              kk = 0
              do j=1,fieldCursor % dimSizes(1)
                k = mpas_binary_search(uniqueSortedNeededList, 2, 1, nUniqueNeededList, fieldCursor % array(j))
                if(k <= nUniqueNeededList) then
                  if (ownerListIn(k) == -i) then
                    kk = kk + 1
                    exchListPtr % srcList(kk) = packingOrder(k)
                    exchListPtr % destList(kk) = j + offsetList(iBlock)
                  end if
                end if
              end do
            end if

            fieldCursor => fieldCursor % next
          end do
        end do

        !
        ! Free up memory
        !
        deallocate(numToSend)
        deallocate(numToRecv)
        deallocate(neededList)
        deallocate(neededListIndex)
        deallocate(neededBlock)

        deallocate(ownedList)
        deallocate(ownedListIndex)
        deallocate(ownedBlock)
        deallocate(ownedListSorted)
        deallocate(ownedBlockSorted)

        deallocate(recipientList)

        deallocate(ownerListIn)
        deallocate(ownerListOut)

        deallocate(uniqueSortedNeededList)
        deallocate(packingOrder)
#endif

        ! Build Copy Lists
        allocate(numToSend(1))
        fieldCursor => ownedListField
        do while (associated(fieldCursor))
          iBlock = fieldCursor % block % localBlockID + 1
          nOwnedList = ownedLimitList(iBlock)
          allocate(ownedListSorted(2, nOwnedList))
          allocate(recipientList(2, nOwnedList))

          do i = 1, nOwnedList
            ownedListSorted(1, i) = fieldCursor % array(i)
            ownedListSorted(2, i) = i
          end do

          call mpas_quicksort(nOwnedList, ownedListSorted)

          fieldCursor2 => neededListField
          do while(associated(fieldCursor2))
            if(associated(fieldCursor, fieldCursor2)) then
              fieldCursor2 => fieldCursor2 % next
              cycle
            end if

            numToSend = 0
            recipientList = -1

            do i = 1, fieldCursor2 % dimSizes(1)
              k = mpas_binary_search(ownedListSorted, 2, 1, nOwnedList, fieldCursor2 % array(i))
              if (k <= nOwnedList) then
                numToSend(1) = numToSend(1) + 1
                ! recipientList(1,:) represents the needed block id
                recipientList(1,ownedListSorted(2,k)) = fieldCursor2 % block % localBlockID
                ! recipientList(2,:) represnets the index in the buffer to place this data
                recipientList(2,ownedListSorted(2,k)) = i
              end if
            end do

            if(numToSend(1) > 0) then
              if(.not.associated(fieldCursor % copyList % halos(haloLayer) % exchList)) then
                allocate(fieldCursor % copyList % halos(haloLayer) % exchList)
                exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList
                nullify(exchListPtr % next)
              else
                ! Find end of copy list
                exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList
                exchListPtr2 => fieldCursor % copyList % halos(haloLayer) % exchList % next
                do while(associated(exchListPtr2))
                  exchListPtr => exchListPtr % next
                  exchListPtr2 => exchListPtr % next
                end do

                allocate(exchListPtr % next)
                exchListPtr => exchListPtr % next
                nullify(exchListPtr % next)
              end if

              exchListPtr % endPointID = fieldCursor2 % block % localBlockID
              exchListPtr % nlist = numToSend(1)
              allocate(exchListPtr % srcList(numToSend(1)))
              allocate(exchListPtr % destList(numToSend(1)))
              exchListPtr % srcList = -1
              exchListPtr % destList = -1

              kk = 1
              do j=1,nOwnedList
               if(recipientList(1,j) == fieldCursor2 % block % localBlockID) then
                 exchListPtr % srcList(kk) = j
                 exchListPtr % destList(kk) = recipientList(2,j) + offSetList(fieldCursor2 % block % localBlockID+1)
                 kk = kk + 1
               end if
              end do
            end if
            fieldCursor2 => fieldCursor2 % next
          end do

          deallocate(recipientList)
          deallocate(ownedListSorted)
          fieldCursor => fieldCursor % next
        end do
        deallocate(numToSend)
        deallocate(offSetList)
        deallocate(ownedLimitList)
      end if

   end subroutine mpas_dmpar_get_exch_list!}}}


!***********************************************************************
!
!  routine mpas_dmpar_build_comm_lists
!
!> \brief   Builds send and receive comm lists templates for populating with buffer data
!> \author  Matt Hoffman
!> \date    8 October 2013
!> \details
!>  This routine builds the templates for send and receive communication lists
!>  that can subsequently be populated with buffer data.  Specifically, it
!>  creates all elements of send and receive linked lists of type
!>  mpas_communication_list and fills in the procID and nlist attributes.
!>  Other dmpar routines can then add data to the rbuffer and/or ibuffer arrays
!>  of these send/receive lists.  This initial step is needed by all of the
!>  various halo-exchange subroutines, so encapsulating these
!>  initial steps here allows significant code-reuse and shortening of those
!>  subroutines.  This implementation avoids the use of 'field' because there are
!>  more than 10 different field types, which prevents the ability to generalize
!>  this routine.  However, since we need to traverse blocks, the subroutine
!>  required adding next/prev pointers to the mpas_multihalo_exchange_list type
!>  (since we can't rely on the field to traverse blocks).  This subroutine
!>  has been made public, so cores have access to it.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_build_comm_lists(sendExchList, recvExchList, haloLayers, dimsizes, sendCommList, recvCommList)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_multihalo_exchange_list), pointer, intent(in) :: &
         sendExchList  !< Input: the send exchange list from the variable for which communication lists are desired
      type (mpas_multihalo_exchange_list), pointer, intent(in) :: &
         recvExchList  !< Input: the receive exchange list from the variable for which communication lists are desired
      integer, dimension(:), pointer :: haloLayers !< Input: list of halo layers to be communicated.
      integer, dimension(:) :: dimSizes !< array of sizes of the dimensions of the variable being communicated

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_communication_list), pointer, intent(inout) :: &
         sendCommList  !< Input/Output: the send communication list, unallocated on input, partially filled out on output
      type (mpas_communication_list), pointer, intent(inout) :: &
         recvCommList  !< Input/Output: the receive communication list, unallocated on input, partially filled out on output

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (mpas_multihalo_exchange_list), pointer :: sendListCursor, recvListCursor
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr, commListPtr2
      logical :: comm_list_found
      integer :: nAdded, bufferOffset
      integer :: iHalo
      integer :: nHaloLayers
      integer :: iDimen
      integer :: dimSizeProduct  ! the product of the size of all dimensions
      integer :: threadNum

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         ! We only need the product of all dimension sizes (except the last), so calculate that now
         dimSizeProduct = 1
         do iDimen = 1, size(dimsizes) - 1
            dimSizeProduct = dimSizeProduct * dimsizes(iDimen)
         enddo

         ! Get size of haloLayers array
         nHaloLayers = size(haloLayers)

         ! Allocate communication lists, and setup dead header node.
         allocate(sendCommList)
         nullify(sendCommList % next)
         sendCommList % procID = -1
         sendCommList % nList = 0

         allocate(recvCommList)
         nullify(recvCommList % next)
         recvCommList % procID = -1
         recvCommList % nList = 0


         ! Determine size of buffers for communication lists
         sendListCursor => sendExchList
         recvListCursor => recvExchList  ! We need to traverse the send and recv exchange lists together in this loop
         do while(associated(sendListCursor))

           ! Need to aggregate across halo layers
           do iHalo = 1, nHaloLayers

             ! Determine size from send lists & build the send list
             exchListPtr => sendListCursor % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))  ! loop through items representing different endPoint Id's
               comm_list_found = .false.

               commListPtr => sendCommList
               do while(associated(commListPtr))  ! Loop through items representing different procs being sent to
                 if(commListPtr % procID == exchListPtr % endPointId) then
                   comm_list_found = .true.
                   commListPtr % nList = commListPtr % nList + exchListPtr % nList * dimSizeProduct
                   exit
                 end if

                 commListPtr => commListPtr % next
               end do

               if(.not. comm_list_found) then  ! Add an item to the sendCommList for this endpoint
                 commListPtr => sendCommList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
                 commListPtr % procID = exchListPtr % endPointID
                 commListPtr % nList = exchListPtr % nList * dimSizeProduct
               end if

               exchListPtr => exchListPtr % next
             end do

             ! Setup recv lists
             exchListPtr => recvListCursor % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               comm_list_found = .false.

               commListPtr => recvCommList
               do while(associated(commListPtr))
                 if(commListPtr % procID == exchListPtr % endPointId) then
                   comm_list_found = .true.
                   commListPtr % nList = commListPtr % nList + exchListPtr % nList * dimSizeProduct
                   exit
                 end if

                 commListPtr => commListPtr % next
               end do

               if(.not. comm_list_found) then
                 commListPtr => recvCommList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
                 commListPtr % procID = exchListPtr % endPointID
                 commListPtr % nList = exchListPtr % nList * dimSizeProduct
               end if

               exchListPtr => exchListPtr % next
             end do
           end do   ! halo loop

           sendListCursor => sendListCursor % next  ! Advance to next block (only happens if more than 1 block per proc)
           recvListCursor => recvListCursor % next  ! Advance to next block (only happens if more than 1 block per proc)
           ! We need to traverse the send and recv exchange lists together in this loop (since we cannot traverse the field itself)
         end do  ! sendListCursor (block loop)

         ! Remove the dead head pointer on send and recv list
         commListPtr => sendCommList
         sendCommList => sendCommList % next
         deallocate(commListPtr)

         commListPtr => recvCommList
         recvCommList => recvCommList % next
         deallocate(commListPtr)

         ! Determine size of receive lists
         commListPtr => recvCommList
         do while(associated(commListPtr))
           bufferOffset = 0
           do iHalo = 1, nHaloLayers
             nAdded = 0

             recvListCursor => recvExchList
             do while(associated(recvListCursor))
               exchListPtr => recvListCursor % halos(haloLayers(iHalo)) % exchList
               do while(associated(exchListPtr))
                 if(exchListPtr % endPointID == commListPtr % procID) then
                   nAdded = max(nAdded, maxval(exchListPtr % srcList) * dimSizeProduct)
                 end if
                 exchListPtr => exchListPtr % next
               end do

               recvListCursor => recvListCursor % next
             end do
             bufferOffset = bufferOffset + nAdded
           end do
           commListPtr % nList = bufferOffset

           commListPtr => commListPtr % next
         end do  ! commListPtr
      end if

   !--------------------------------------------------------------------
   end subroutine mpas_dmpar_build_comm_lists!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field1d_integer
!
!> \brief MPAS dmpar all-to-all 1D integer routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field1d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field1dInteger), pointer :: fieldIn !< Input: Field to send
     type (field1dInteger), pointer :: fieldOut !< Output: Field to receive
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: Halo layers to communicated. Defaults to all.

     type (field1dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     threadNum = mpas_threading_get_thread_num()
     dminfo => fieldIn % block % domain % dminfo

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = 0
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         commListPtr % ibuffer = 0
         call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   iBuffer = exchListPtr % destList(i) + bufferOffset
                   commListPtr % ibuffer(iBuffer) = fieldInPtr % array(exchListPtr % srcList(i))
                   nAdded = nAdded + 1
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(exchListPtr % destList(i)) = fieldInPtr % array(exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   iBuffer = exchListPtr % srcList(i) + bufferOffset
                   fieldOutPtr % array(exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field1d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field2d_integer
!
!> \brief MPAS dmpar all-to-all 2D integer routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field2d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field2dInteger), pointer :: fieldIn !< Input: Field to communicate from
     type (field2dInteger), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field2dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = bufferOffset

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(1)
                     iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) + j + bufferOffset
                     commListPtr % ibuffer(iBuffer) = fieldInPtr % array(j, exchListPtr % srcList(i))
                     nAdded = nAdded + 1
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, exchListPtr % destList(i)) = fieldInPtr % array(:, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(1)
                     iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(1) + j + bufferOffset
                     fieldOutPtr % array(j, exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field2d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field3d_integer
!
!> \brief MPAS dmpar all-to-all 3D integer routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field3d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field3dInteger), pointer :: fieldIn !< Input: Field to send from
     type (field3dInteger), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field3dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers
       commListPtr => recvList
       do while(associated(commListPtr))
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % ibuffer(commListPtr % nList))
         nullify(commListPtr % rbuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(2)
                     do k = 1, fieldInPtr % dimSizes(1)
                       iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) + (j-1) * fieldInPtr % dimSizes(1) + k + bufferOffset
                       commListPtr % ibuffer(iBuffer) = fieldInPtr % array(k, j, exchListPtr % srcList(i))
                       nAdded = nAdded + 1
                     end do
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(2)
                     do k = 1, fieldOutPtr % dimSizes(1)
                       iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1) + (j-1) * fieldOutPtr % dimSizes(1) + k + bufferOffset
                       fieldOutPtr % array(k, j, exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
                     end do
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field3d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field1d_real
!
!> \brief MPAS dmpar all-to-all 1D real routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field1d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field1dReal), pointer :: fieldIn !< Input: Field to send from
     type (field1dReal), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field1dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers
       commListPtr => recvList
       do while(associated(commListPtr))
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   iBuffer = exchListPtr % destList(i) + bufferOffset
                   commListPtr % rbuffer(iBuffer) = fieldInPtr % array(exchListPtr % srcList(i))
                   nAdded = nAdded + 1
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(exchListPtr % destList(i)) = fieldInPtr % array(exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   iBuffer = exchListPtr % srcList(i) + bufferOffset
                   fieldOutPtr % array(exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field1d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field2d_real
!
!> \brief MPAS dmpar all-to-all 2D real routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field2d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field2dReal), pointer :: fieldIn !< Input: Field to send from
     type (field2dReal), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field2dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(1)
                     iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) + j + bufferOffset
                     commListPtr % rbuffer(iBuffer) = fieldInPtr % array(j, exchListPtr % srcList(i))
                     nAdded = nAdded + 1
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, exchListPtr % destList(i)) = fieldInPtr % array(:, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(1)
                     iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(1) + j + bufferOffset
                     fieldOutPtr % array(j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field2d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field3d_real
!
!> \brief MPAS dmpar all-to-all 3D real routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field3d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field3dReal), pointer :: fieldIn !< Input: Field to send from
     type (field3dReal), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field3dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(2)
                     do k = 1, fieldInPtr % dimSizes(1)
                       iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) + (j-1) * fieldInPtr % dimSizes(1) + k + bufferOffset
                       commListPtr % rbuffer(iBuffer) = fieldInPtr % array(k, j, exchListPtr % srcList(i))
                       nAdded = nAdded + 1
                     end do
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(2)
                     do k = 1, fieldOutPtr % dimSizes(1)
                       iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1) + (j-1) * fieldOutPtr % dimSizes(1) + k + bufferOffset
                       fieldOutPtr % array(k, j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                     end do
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field3d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field4d_real
!
!> \brief MPAS dmpar all-to-all 4D real routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field4d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field4dReal), pointer :: fieldIn !< Input: Field to send from
     type (field4dReal), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

     type (field4dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k, l
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldOutPtr % dimSizes(3)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldOutPtr % dimSizes(3)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(3)
                     do k = 1, fieldInPtr % dimSizes(2)
                       do l = 1, fieldInPtr % dimSizes(1)
                         iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldInPtr % dimSizes(3) &
                                 + (j-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) &
                                 + (k-1) * fieldInPtr % dimSizes(1) + l + bufferOffset
                         commListPtr % rbuffer(iBuffer) = fieldInPtr % array(l, k, j, exchListPtr % srcList(i))
                         nAdded = nAdded + 1
                       end do
                     end do
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, :, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, :, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(3)
                     do k = 1, fieldOutPtr % dimSizes(2)
                       do l = 1, fieldOutPtr % dimSizes(1)
                         iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1)  &
                                 + (j-1) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)  &
                                 + (k-1) * fieldOutPtr % dimSizes(1) + l + bufferOffset
                         fieldOutPtr % array(l, k, j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                       end do
                     end do
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field4d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_alltoall_field5d_real
!
!> \brief MPAS dmpar all-to-all 5D real routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the all-to-all communication of an input field into an output field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_alltoall_field5d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field5dReal), pointer :: fieldIn !< Input: Field to send from
     type (field5dReal), pointer :: fieldOut !< Output: Field to receive into
     integer, dimension(:), pointer, optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all.

     type (field5dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k, l, m
     integer :: nHaloLayers, threadNum
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       if(present(haloLayersIn)) then
         nHaloLayers = size(haloLayersIn)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = haloLayersIn(iHalo)
         end do
       else
         nHaloLayers = size(fieldIn % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
           haloLayers(iHalo) = iHalo
         end do
       end if

#ifdef _MPI
       nullify(sendList)
       nullify(recvList)

       ! Setup receive lists.
       do iHalo = 1, nHaloLayers
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => recvList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(4)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(recvList)) then
                 allocate(recvList)
                 nullify(recvList % next)
                 commListPtr => recvList
               else
                 commListPtr => recvList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if

               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(4)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldOutPtr => fieldOutPtr % next
         end do
       end do

       ! Determine size of receive list buffers.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(4))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do
         commListPtr % nList = nAdded

         commListPtr => commListPtr % next
       end do

       ! Allocate buffers for receives, and initiate mpi_irecv calls.
       commListPtr => recvList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Setup send lists, and determine the size of their buffers.
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             comm_list_found = .false.

             ! Search for an already created commList to this processor.
             commListPtr => sendList
             do while(associated(commListPtr))
               if(commListPtr % procID == exchListPtr % endPointID) then
                 commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldInPtr % dimSizes(4)
                 comm_list_found = .true.
                 exit
               end if

               commListPtr => commListPtr % next
             end do

             ! If no comm list exists, create a new one.
             if(.not. comm_list_found) then
               if(.not.associated(sendList)) then
                 allocate(sendList)
                 nullify(sendList % next)
                 commListPtr => sendList
               else
                 commListPtr => sendList
                 commListPtr2 => commListPtr % next
                 do while(associated(commListPtr2))
                   commListPtr => commListPtr % next
                   commListPtr2 => commListPtr % next
                 end do

                 allocate(commListPtr % next)
                 commListPtr => commListPtr % next
                 nullify(commListPtr % next)
               end if
               commListPtr % procID = exchListPtr % endPointID
               commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldInPtr % dimSizes(4)
             end if

             exchListPtr => exchListPtr % next
           end do

           fieldInPtr => fieldInPtr % next
         end do
       end do

       ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
       commListPtr => sendList
       do while(associated(commListPtr))
         allocate(commListPtr % rbuffer(commListPtr % nList))
         nullify(commListPtr % ibuffer)
         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldInPtr => fieldIn
           do while(associated(fieldInPtr))
             exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldInPtr % dimSizes(4)
                     do k = 1, fieldInPtr % dimSizes(3)
                       do l = 1, fieldInPtr % dimSizes(2)
                         do m = 1, fieldInPtr % dimSizes(1)
                           iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldInPtr % dimSizes(3) * fieldInPtr % dimSizes(4) &
                                   + (j-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) * fieldInPtr % dimSizes(3) &
                                   + (k-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) &
                                   + (l-1) * fieldInPtr % dimSizes(1) + m + bufferOffset
                           commListPtr % rbuffer(iBuffer) = fieldInPtr % array(m, l, k, j, exchListPtr % srcList(i))
                           nAdded = nAdded + 1
                         end do
                       end do
                     end do
                   end do
                 end do
               end if

               exchListPtr => exchListPtr % next
             end do

             fieldInPtr => fieldInPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                        commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

         commListPtr => commListPtr % next
       end do

#endif

       ! Handle Local Copies. Only local copies if no MPI
       do iHalo = 1, nHaloLayers
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             fieldOutPtr => fieldOut
             do while(associated(fieldOutPtr))
               if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
                 do i = 1, exchListPtr % nList
                   fieldOutPtr % array(:, :, :, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, :, :, exchListPtr % srcList(i))
                 end do
               end if
               fieldOutPtr => fieldOutPtr % next
             end do

             exchListPtr => exchListPtr % next
           end do
           fieldInPtr => fieldInPtr % next
         end do
       end do

#ifdef _MPI
       ! Wait for MPI_Irecv's to finish, and unpack data.
       commListPtr => recvList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

         bufferOffset = 0
         do iHalo = 1, nHaloLayers
           nAdded = 0
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
             do while(associated(exchListPtr))
               if(exchListPtr % endPointID == commListPtr % procID) then
                 do i = 1, exchListPtr % nList
                   do j = 1, fieldOutPtr % dimSizes(4)
                     do k = 1, fieldOutPtr % dimSizes(3)
                       do l = 1, fieldOutPtr % dimSizes(2)
                         do m = 1, fieldOutPtr % dimSizes(1)
                           iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(4) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1) &
                                   + (j-1) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) &
                                   + (k-1) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) &
                                   + (l-1) * fieldOutPtr % dimSizes(1) + m + bufferOffset
                           fieldOutPtr % array(m, l, k, j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                         end do
                       end do
                     end do
                   end do
                 end do
                 nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(3) * fieldOutPtr % dimSizes(4))
               end if
               exchListPtr => exchListPtr % next
             end do

             fieldOutPtr => fieldOutPtr % next
           end do
           bufferOffset = bufferOffset + nAdded
         end do

         commListPtr => commListPtr % next
       end do

       ! Wait for MPI_Isend's to finish.
       commListPtr => sendList
       do while(associated(commListPtr))
         call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         commListPtr => commListPtr % next
       end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_alltoall_field5d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field1d_integer
!
!> \brief MPAS dmpar halo exchange 1D integer field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field1d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field1DInteger), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field1DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    commListPtr % ibuffer(exchListPtr % destList(i) + bufferOffset) = fieldCursor % array(exchListPtr % srcList(i))
                    nAdded = nAdded + 1

                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(exchListPtr % destList(i)) = fieldCursor % array(exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor % array(exchListPtr % destList(i)) = commListPtr % ibuffer(exchListPtr % srcList(i) + bufferOffset)
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

        ! Destroy commLists.
        call mpas_dmpar_destroy_communication_list(sendList)
        call mpas_dmpar_destroy_communication_list(recvList)
#endif

        deallocate(haloLayers)
      end if

   end subroutine mpas_dmpar_exch_halo_field1d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field2d_integer
!
!> \brief MPAS dmpar halo exchange 2D integer field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field2d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field2DInteger), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field2DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(1)
                      commListPtr % ibuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) + j + bufferOffset) = fieldCursor % array(j, exchListPtr % srcList(i))
                      nAdded = nAdded + 1
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, exchListPtr % destList(i)) = fieldCursor % array(:, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(1)
                      fieldCursor % array(j, exchListPtr % destList(i)) = commListPtr % ibuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) + j + bufferOffset)
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field2d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field3d_integer
!
!> \brief MPAS dmpar halo exchange 3D integer field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field3d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field3DInteger), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field3DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadnum
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % ibuffer(commListPtr % nList))
          nullify(commListPtr % rbuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(2)
                      do k = 1, fieldCursor % dimSizes(1)
                        commListPtr % ibuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                            + (j-1) * fieldCursor % dimSizes(1) + k  + bufferOffset) = fieldCursor % array(k, j, exchListPtr % srcList(i))
                        nAdded = nAdded + 1
                      end do
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(2)
                      do k = 1, fieldCursor % dimSizes(1)
                        fieldCursor % array(k, j, exchListPtr % destList(i)) = commListPtr % ibuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                             + (j-1)*fieldCursor % dimSizes(1) + k + bufferOffset)
                      end do
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field3d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field1d_real
!
!> \brief MPAS dmpar halo exchange 1D real field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field1d_real(field, haloLayersIn)!{{{

      implicit none

      type (field1dReal), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field1dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    commListPtr % rbuffer(exchListPtr % destList(i) + bufferOffset) = fieldCursor % array(exchListPtr % srcList(i))
                    nAdded = nAdded + 1
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(exchListPtr % destList(i)) = fieldCursor % array(exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor % array(exchListPtr % destList(i)) = commListPtr % rbuffer(exchListPtr % srcList(i) + bufferOffset)
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field1d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field2d_real
!
!> \brief MPAS dmpar halo exchange 2D real field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field2d_real(field, haloLayersIn)!{{{

      implicit none

      type (field2dReal), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field2dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI
        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(1)
                      commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) + j + bufferOffset) = fieldCursor % array(j, exchListPtr % srcList(i))
                      nAdded = nAdded + 1
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, exchListPtr % destList(i)) = fieldCursor % array(:, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(1)
                      fieldCursor % array(j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizeS(1) + j + bufferOffset)
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field2d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field3d_real
!
!> \brief MPAS dmpar halo exchange 3D real field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field3d_real(field, haloLayersIn)!{{{

      implicit none

      type (field3dReal), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field3dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI
        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(2)
                      do k = 1, fieldCursor % dimSizes(1)
                        commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                            + (j-1) * fieldCursor % dimSizes(1) + k  + bufferOffset) = fieldCursor % array(k, j, exchListPtr % srcList(i))
                        nAdded = nAdded + 1
                      end do
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(2)
                      do k = 1, fieldCursor % dimSizes(1)
                        fieldCursor % array(k, j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                             + (j-1)*fieldCursor % dimSizes(1) + k + bufferOffset)
                      end do
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field3d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field4d_real
!
!> \brief MPAS dmpar halo exchange 4D real field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field4d_real(field, haloLayersIn)!{{{

      implicit none

      type (field4dReal), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field4dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i, j, k, l
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 4
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(3)
                      do k = 1, fieldCursor % dimSizes(2)
                        do l = 1, fieldCursor % dimSizes(1)
                          commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                              + (j-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                              + (k-1) * fieldCursor % dimSizes(1) + l  + bufferOffset) &
                              = fieldCursor % array(l, k, j, exchListPtr % srcList(i))
                          nAdded = nAdded + 1
                        end do
                      end do
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, :, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, :, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(3)
                      do k = 1, fieldCursor % dimSizes(2)
                        do l = 1, fieldCursor % dimSizes(1)
                          fieldCursor % array(l, k, j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)&
                                                                                 *fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) *fieldCursor % dimSizes(3)&
                                                                               + (j-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                               + (k-1)*fieldCursor % dimSizes(1) + l + bufferOffset)
                        end do
                      end do
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field4d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_halo_field5d_real
!
!> \brief MPAS dmpar halo exchange 5D real field
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine handles the halo exchange communication of an input field across all processors.
!>  It requires exchange lists to be created prior to calling this routine.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_halo_field5d_real(field, haloLayersIn)!{{{

      implicit none

      type (field5dReal), pointer :: field !< Input: Field to communicate
      integer, dimension(:), intent(in), optional :: haloLayersIn !< Input: List of halo layers to communicate. Defaults to all

      type (dm_info), pointer :: dminfo
      type (field5dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr
      integer :: mpi_ierr, threadNum
      integer :: nHaloLayers, iHalo, i, j, k, l, m
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      if ( .not. field % isActive ) then
#ifdef MPAS_DEBUG
         write(stderrUnit, *) ' -- Skipping halo exchange for deactivated field: ' // trim(field % fieldName)
#endif
         return
      end if

      do i = 1, 5
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo
      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
        if(present(haloLayersIn)) then
          nHaloLayers = size(haloLayersIn)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = haloLayersIn(iHalo)
          end do
        else
          nHaloLayers = size(field % sendList % halos)
          allocate(haloLayers(nHaloLayers))
          do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
          end do
        end if

#ifdef _MPI

        ! Setup Communication Lists
        call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendList, recvList)

        ! Allocate space in recv lists, and initiate mpi_irecv calls
        commListPtr => recvList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

          commListPtr => commListPtr % next
        end do

        ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
        commListPtr => sendList
        do while(associated(commListPtr))
          allocate(commListPtr % rbuffer(commListPtr % nList))
          nullify(commListPtr % ibuffer)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do  i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(4)
                      do k = 1, fieldCursor % dimSizes(3)
                        do l = 1, fieldCursor % dimSizes(2)
                          do m = 1, fieldCursor % dimSizes(1)
                            commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4) &
                                + (j-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                + (k-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                + (l-1) * fieldCursor % dimSizes(1) + m + bufferOffset) &
                                = fieldCursor % array(m, l, k, j, exchListPtr % srcList(i))
                            nAdded = nAdded + 1
                          end do
                        end do
                      end do
                    end do
                  end do
                end if

                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do

          call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
          commListPtr => commListPtr % next
        end do
#endif

        ! Handle local copy. If MPI is off, then only local copies are performed.
        fieldCursor => field
        do while(associated(fieldCursor))
          do iHalo = 1, nHaloLayers
            exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

            do while(associated(exchListPtr))
              fieldCursor2 => field
              do while(associated(fieldCursor2))
                if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  do i = 1, exchListPtr % nList
                    fieldCursor2 % array(:, :, :, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, :, :, exchListPtr % srcList(i))
                  end do
                end if

                fieldCursor2 => fieldCursor2 % next
              end do

              exchListPtr => exchListPtr % next
            end do
          end do

          fieldCursor => fieldCursor % next
        end do

#ifdef _MPI

        ! Wait for mpi_irecv to finish, and unpack data from buffer
        commListPtr => recvList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          bufferOffset = 0
          do iHalo = 1, nHaloLayers
            nAdded = 0
            fieldCursor => field
            do while(associated(fieldCursor))
              exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
              do while(associated(exchListPtr))
                if(exchListPtr % endPointID == commListPtr % procID) then
                  do i = 1, exchListPtr % nList
                    do j = 1, fieldCursor % dimSizes(4)
                      do k = 1, fieldCursor % dimSizes(3)
                        do l = 1, fieldCursor % dimSizes(2)
                          do m = 1, fieldCursor % dimSizes(1)
                            fieldCursor % array(m, l, k, j, exchListPtr % destList(i)) = &
                                commListPtr % rbuffer((exchListPtr % srcList(i)-1) &
                                  * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                    * fieldCursor % dimSizes(4)&
                                  + (j-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                  + (k-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                  + (l-1)*fieldCursor % dimSizes(1) + m + bufferOffset)
                          end do
                        end do
                      end do
                    end do
                  end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                               * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4))
                end if
                exchListPtr => exchListPtr % next
              end do

              fieldCursor => fieldCursor % next
            end do
            bufferOffset = bufferOffset + nAdded
          end do
          commListPtr => commListPtr % next
        end do

        ! wait for mpi_isend to finish.
        commListPtr => sendList
        do while(associated(commListPtr))
          call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
          commListPtr => commListPtr % next
        end do

       ! Destroy commLists.
       call mpas_dmpar_destroy_communication_list(sendList)
       call mpas_dmpar_destroy_communication_list(recvList)
#endif

       deallocate(haloLayers)
     end if

   end subroutine mpas_dmpar_exch_halo_field5d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_init_multihalo_exchange_list
!
!> \brief MPAS dmpar initialize muiltihalo exchange list routine.
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine initializes the multihalo exchange lists, based on a number of halo layers.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_init_multihalo_exchange_list(exchList, nHalos)!{{{
     type (mpas_multihalo_exchange_list), pointer :: exchList !< Input: Exchange list to initialize
     integer, intent(in) :: nHalos !< Input: Number of halo layers for exchange list

     integer :: i, threadNum

     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       allocate(exchList)
       allocate(exchList % halos(nHalos))
       do i = 1, nHalos
         nullify(exchList % halos(i) % exchList)
       end do
       nullify(exchList % next)
       nullify(exchList % prev)
     end if
   end subroutine mpas_dmpar_init_multihalo_exchange_list!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_destroy_mulithalo_exchange_list
!
!> \brief MPAS dmpar destroy muiltihalo exchange list routine.
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine destroys the multihalo exchange lists.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_destroy_mulithalo_exchange_list(exchList)!{{{
     type (mpas_multihalo_exchange_list), pointer :: exchList !< Input: Exchange list to destroy.

     integer :: nHalos
     integer :: i, threadNum

     nHalos = size(exchList % halos)
     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       do i = 1, nHalos
         call mpas_dmpar_destroy_exchange_list(exchList % halos(i) % exchList)
       end do

       deallocate(exchList % halos)
       deallocate(exchList)
       nullify(exchList)
     end if
   end subroutine mpas_dmpar_destroy_mulithalo_exchange_list!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_destroy_communication_list
!
!> \brief MPAS dmpar destroy communication list routine.
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine destroys a communication lists.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_destroy_communication_list(commList)!{{{
     type (mpas_communication_list), pointer :: commList !< Input: Communication list to destroy.
     type (mpas_communication_list), pointer :: commListPtr
     integer :: threadNum

     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       commListPtr => commList
       do while(associated(commListPtr))
         if(associated(commList % next)) then
           commList => commList % next
         else
           nullify(commList)
         end if

         if(associated(commListPtr % ibuffer)) then
           deallocate(commListPtr % ibuffer)
         end if

         if(associated(commListPtr % rbuffer)) then
           deallocate(commListPtr % rbuffer)
         end if

         deallocate(commListPtr)
         commListPtr => commList
       end do
     end if

   end subroutine mpas_dmpar_destroy_communication_list!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_destroy_exchange_list
!
!> \brief MPAS dmpar destroy exchange list routine.
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine destroys a exchange lists.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_destroy_exchange_list(exchList)!{{{
     type (mpas_exchange_list), pointer :: exchList !< Input: Exchange list to destroy
     type (mpas_exchange_list), pointer :: exchListPtr
     integer :: threadNum

     threadNum = mpas_threading_get_thread_num()

     if ( threadNum == 0 ) then
       exchListPtr => exchList
       do while(associated(exchList))
         if(associated(exchList % next)) then
           exchList => exchList % next
         else
           nullify(exchList)
         end if

         if(associated(exchListPtr % srcList)) then
           deallocate(exchListPtr % srcList)
         end if

         if(associated(exchListPtr % destList)) then
           deallocate(exchListPtr % destList)
         end if

         deallocate(exchListPtr)
         exchListPtr => exchList
       end do
     end if

   end subroutine mpas_dmpar_destroy_exchange_list!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field1d_integer
!
!> \brief MPAS dmpar copy 1D integer field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 1D integer field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field1d_integer(field)!{{{
       type (field1dInteger), pointer :: field !< Input: Field to copy
       type (field1dInteger), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field1d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field2d_integer
!
!> \brief MPAS dmpar copy 2D integer field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 2D integer field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field2d_integer(field)!{{{
       type (field2dInteger), pointer :: field !< Input: Field to copy
       type (field2dInteger), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field2d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field3d_integer
!
!> \brief MPAS dmpar copy 3D integer field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 3D integer field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field3d_integer(field)!{{{
       type (field3dInteger), pointer :: field !< Input: Field to copy
       type (field3dInteger), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field3d_integer!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field1d_real
!
!> \brief MPAS dmpar copy 1D real field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 1D real field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field1d_real(field)!{{{
       type (field1dReal), pointer :: field !< Input: Field to copy
       type (field1dReal), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field
           do while(associated(fieldCursor))
             fieldCursor % array(:) = field % array(:)
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field1d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field2d_real
!
!> \brief MPAS dmpar copy 2D real field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 2D real field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field2d_real(field)!{{{
       type (field2dReal), pointer :: field !< Input: Field to copy
       type (field2dReal), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field2d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field3d_real
!
!> \brief MPAS dmpar copy 3D real field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 3D real field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field3d_real(field)!{{{
       type (field3dReal), pointer :: field !< Input: Field to copy
       type (field3dReal), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field3d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field4d_real
!
!> \brief MPAS dmpar copy 4D real field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 4D real field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field4d_real(field)!{{{
       type (field4dReal), pointer :: field !< Input: Field to copy
       type (field4dReal), pointer :: fieldCursor
       integer :: threadNum
       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field4d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_copy_field5d_real
!
!> \brief MPAS dmpar copy 5D real field routine
!> \author Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine copies a 5D real field throughout a block list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_copy_field5d_real(field)!{{{
       type (field5dReal), pointer :: field !< Input: Field to copy
       type (field5dReal), pointer :: fieldCursor
       integer :: threadNum

       threadNum = mpas_threading_get_thread_num()

       if ( threadNum == 0 ) then
         if(associated(field % next)) then
           fieldCursor => field % next
           do while(associated(fieldCursor))
             fieldCursor % array = field % array
             fieldCursor => fieldCursor % next
           end do
         end if
       end if
   end subroutine mpas_dmpar_copy_field5d_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_get_time
!
!> \brief MPAS dmpar get time routine
!> \author Doug Jacobsen
!> \date   06/11/2015
!> \details
!>  This routine returns the current time, either using the MPI interface, or
!>  the system_clock interface.
!
!-----------------------------------------------------------------------

   subroutine mpas_dmpar_get_time(curTime)!{{{

      implicit none

      real (kind=R8KIND), intent(out) :: curTime !< Output: Current time
      integer :: clock, hz, threadNum

#ifdef _MPI
      curTime = MPI_WTime()
#else
      call system_clock(count=clock)
      call system_clock(count_rate=hz)
      curTime = real(clock, kind=R8KIND) / real(hz, kind=R8KIND)
#endif

   end subroutine mpas_dmpar_get_time!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_create
!
!> \brief MPAS dmpar exchange group creation routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine creates a new exchange group named 'groupName' within a
!>  domain's list of exchange groups.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_create(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_group), pointer :: exchGroupPtr, prevExchGroupPtr
      integer :: nLen
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         nLen = len_trim(groupName)
         DMPAR_DEBUG_WRITE(' -- Creating exchange group ' // trim(groupName))

         ! If no exchange group has been created yet, allocate the first
         if ( .not. associated(domain % exchangeGroups) ) then
            allocate(domain % exchangeGroups)
            exchGroupPtr => domain % exchangeGroups
         ! Otherwise, make sure the group has not already been created, and add it to the end of the list
         else
            exchGroupPtr => domain % exchangeGroups
            do while (associated(exchGroupPtr))
               if ( nLen == exchGroupPtr % nLen) then
                  if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
                     write(stderrUnit, *) 'ERROR: Exchange group ''' // groupName(1:nLen) // ''' already exists.'
                     write(stderrUnit, *) '       Cannot create group.'
                     if ( present(iErr) ) then
                        iErr = MPAS_DMPAR_EXISTING_GROUP
                     end if
                     return
                  end if
               end if
               prevExchGroupPtr => exchGroupPtr
               exchGroupPtr => exchGroupPtr % next
            end do

            allocate(prevExchGroupPtr % next)
            exchGroupPtr => prevExchGroupPtr % next
         end if

         exchGroupPtr % nLen = nLen
         exchGroupPtr % groupName = trim(groupName)

         call mpas_pool_create_pool(exchGroupPtr % fieldPool)
      end if

   end subroutine mpas_dmpar_exch_group_create!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_destroy
!
!> \brief MPAS dmpar exchange group destruction routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine destroys an exchange group named 'groupName' within a
!>  domain's list of exchange groups.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_destroy(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr
      type (mpas_exchange_group), pointer :: exchGroupPtr, prevExchGroupPtr
      integer :: nLen
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         nLen = len_trim(groupName)
         DMPAR_DEBUG_WRITE(' -- Destroying exchange group ' // trim(groupName))

         nullify(prevExchGroupPtr)
         exchGroupPtr => domain % exchangeGroups
         do while (associated(exchGroupPtr))
            if ( nLen == exchGroupPtr % nLen) then
               if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
                  if ( associated(prevExchGroupPtr) ) then
                     prevExchGroupPtr % next => exchGroupPtr % next
                  else
                     if ( associated(exchGroupPtr % next) ) then
                        domain % exchangeGroups => exchGroupPtr % next
                     else
                        nullify(domain % exchangeGroups)
                     end if
                  end if

                  call mpas_pool_destroy_pool(exchGroupPtr % fieldPool)

                  ! Destroy field list
                  do while ( associated(exchGroupPtr % fieldList) )
                     exchFieldListPtr => exchGroupPtr % fieldList
                     if ( associated(exchFieldListPtr % next) ) then
                        exchGroupPtr % fieldList => exchFieldListPtr % next
                     else
                        nullify(exchGroupPtr % fieldList)
                     end if

                     if ( associated(exchFieldListPtr % haloLayers) ) then
                        deallocate(exchFieldListPtr % haloLayers)
                     end if

                     if ( associated(exchFieldListPtr % timeLevels) ) then
                        deallocate(exchFieldListPtr % timeLevels)
                     end if

                     nullify(exchFieldListPtr % next)
                     deallocate(exchFieldListPtr)
                  end do

                  deallocate(exchGroupPtr)

                  return
               end if
            end if
            prevExchGroupPtr => exchGroupPtr
            exchGroupPtr => exchGroupPtr % next
         end do

         write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
         write(stderrUnit, *) '       Cannot destroy group.'
         if ( present(iErr) ) then
            iErr = MPAS_DMPAR_MISSING_GROUP
         end if
      end if

   end subroutine mpas_dmpar_exch_group_destroy!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_add_field
!
!> \brief MPAS dmpar exchange group add field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine adds a field named 'fieldName' to an exchange group named
!>  'groupName' within a domain's list of exchange groups.
!>  The timeLevel argument allows control over which timeLevel(s) will be
!>    exchanged as part of this group. If the timeLevel argument is omitted or if it
!>    has a value of -1, all time levels will be exchanged.
!>  The haloLayers argument allows an input array to define the halo layers
!>    that should be exchanged as part of this exchange group. If it is
!>    omitted, all halo layers will be exchanged.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_add_field(domain, groupName, fieldName, timeLevel, haloLayers, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      character (len=*), intent(in) :: fieldName
      integer, optional, intent(in) :: timeLevel
      integer, dimension(:), optional, intent(in) :: haloLayers
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr
      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen, fieldNLen, timeLevelLocal
      type (mpas_pool_field_info_type) :: fieldInfo
      integer :: iHalo
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         nLen = len_trim(groupName)
         fieldNLen = len_trim(fieldName)

         DMPAR_DEBUG_WRITE(' -- Adding field ' // trim(fieldName) // ' to exchange group ' // trim(groupName) )

         ! Check that field exists
         fieldInfo % fieldType = -1
         fieldInfo % nDims = -1
         fieldInfo % nTimeLevels = -1
         fieldInfo % isActive = .false.
         call mpas_pool_get_field_info(domain % blocklist % allFields, fieldName, fieldInfo)

         if ( fieldInfo % fieldType == -1 .and. fieldInfo % nDims == -1 .and. fieldInfo % nTimeLevels == -1 ) then
            write(stderrUnit, *) 'ERROR: Field ''' // trim(fieldName) // ''' does not exist.'
            write(stderrUnit, *) '       Cannot add field.'
            if ( present(iErr) ) then
               iErr = MPAS_DMPAR_MISSING_FIELD
            end if
            return
         end if

         ! Validate dimensions
         if ( fieldInfo % nDims == 0 ) then
            write(stderrUnit, *) 'WARNING: Field ''' // trim(fieldName) // &
                                 ''' has zero dimensions and cannot be added to exchange group ''' // &
                                 trim(groupName) // '''.'
            return
         end if

         ! Validate field type
         if ( fieldInfo % fieldType /= MPAS_POOL_REAL .and. fieldInfo % fieldType /= MPAS_POOL_INTEGER ) then
            write(stderrUnit, *) 'WARNING: Field ''' // trim(fieldName) // &
                                 ''' is an incorrect type for a halo exchange, and cannot be added to exchange group ''' // &
                                 trim(groupName) // '''.'
            return
         end if

         ! Setup and validate time level
         if ( present(timeLevel) ) then
            timeLevelLocal = timeLevel
         else
            timeLevelLocal = -1
         end if

         if ( timeLevelLocal > fieldInfo % nTimeLevels ) then
            write(stderrUnit, *) 'ERROR: Field ''' // trim(fieldName) // &
                                 ''' has fewer time levels than the level requested. Cannot add to exchange group ''' // &
                                trim(groupName) // '''.'
            if ( present(iErr) ) then
               iErr = MPAS_DMPAR_FIELD_TIMELEVEL_ERR
            end if
            return
         else if ( timeLevelLocal < -1 ) then
            write(stderrUnit, '(a62, i2, a)') 'ERROR: timeLevel argument can only have a value between -1 and ', &
                                             fieldInfo % nTimeLevels, ' when adding field ''' // trim(fieldName) &
                                             // ''' to exchange group ''' // trim(groupName) // '''.'
            if ( present(iErr) ) then
               iErr = MPAS_DMPAR_FIELD_TIMELEVEL_ERR
            end if
            return
         end if

         if ( present(haloLayers) ) then
            do iHalo = 1, size(haloLayers)
               if ( haloLayers(iHalo) > fieldInfo % nHaloLayers ) then
                  write(stderrUnit, '(a60, i2, a)') 'ERROR: haloLayers argument contains an invalid halo index of ', &
                                                    haloLayers(iHalo), ' when adding field ''' // trim(fieldName) // &
                                                    ''' to exchange group ''' // trim(groupName) // '''.'
                  if ( present(iErr) ) then
                     iErr = MPAS_DMPAR_FIELD_HALO_ERR
                  end if
                  return
               end if
            end do
         end if

         exchGroupPtr => domain % exchangeGroups
         do while (associated(exchGroupPtr))
            if ( nLen == exchGroupPtr % nLen) then
               if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
                  call mpas_pool_add_config(exchGroupPtr % fieldPool, fieldName, timeLevelLocal)

                  ! Add field to exchange group's fieldList
                  if ( associated(exchGroupPtr % fieldList) ) then
                     ! If the field is already added, find it.
                     exchFieldListPtr => exchGroupPtr % fieldList
                     do while ( associated(exchFieldListPtr) )
                        if ( fieldNLen == exchFieldListPtr % nLen )  then
                           if ( fieldName(1:fieldNLen) == exchFieldListPtr % fieldName(1:exchFieldListPtr % nLen) ) then
                              ! Make sure this field is marked for the correct
                              ! halo exchanges, as defined by the input
                              ! arguments
                              ! TODO
                              if ( timeLevelLocal == -1 ) then
                                 exchFieldListPtr % timeLevels(:) = .true.
                              else
                                 exchFieldListPtr % timeLevels(timeLevel) = .true.
                              end if

                              if ( present(haloLayers) ) then
                                 do iHalo = 1, size(haloLayers)
                                    exchFieldListPtr % haloLayers( haloLayers(iHalo) ) = .true.
                                 end do
                              else
                                 exchFieldListPtr % haloLayers(:) = .true.
                              end if

                              ! Return, as the field has already been added
                              return
                           end if
                        end if

                        exchFieldListPtr => exchFieldListPtr % next
                     end do
                  end if

                  ! Add the field to the beginning of the list, since it
                  ! isn't in the list.
                  allocate(exchFieldListPtr)
                  exchFieldListPtr % nLen = fieldNLen
                  exchFieldListPtr % fieldName = trim(fieldName)
                  exchFieldListPtr % nDims = fieldInfo % nDims
                  exchFieldListPtr % fieldType = fieldInfo % fieldType
                  allocate(exchFieldListPtr % haloLayers(fieldInfo % nHaloLayers))
                  allocate(exchFieldListPtr % timeLevels(fieldInfo % nTimeLevels))
                  exchFieldListPtr % haloLayers(:) = .false.
                  exchFieldListPtr % timeLevels(:) = .false.

                  if ( timeLevelLocal == -1 ) then
                     exchFieldListPtr % timeLevels(:) = .true.
                  else
                     exchFieldListPtr % timeLevels(timeLevelLocal) = .true.
                  end if

                  if ( present(haloLayers) ) then
                     do iHalo = 1, size(haloLayers)
                        exchFieldListPtr % haloLayers( haloLayers(iHalo) ) = .true.
                     end do
                  else
                     exchFieldListPtr % haloLayers(:) = .true.
                  end if

                  ! The next pointer is null by default, so only update it if
                  ! there is a fieldList already.
                  if ( associated(exchGroupPtr % fieldList) ) then
                     exchFieldListPtr % next => exchGroupPtr % fieldList
                  end if
                  exchGroupPtr % fieldList => exchFieldListPtr

                  return
               end if
            end if
            exchGroupPtr => exchGroupPtr % next
         end do

         write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
         write(stderrUnit, *) '       Cannot add field ''' // trim(fieldName) // ''' to group.'
         if ( present(iErr) ) then
            iErr = MPAS_DMPAR_MISSING_GROUP
         end if
      end if

   end subroutine mpas_dmpar_exch_group_add_field!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_remove_field
!
!> \brief MPAS dmpar exchange group remove field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine removes a field named 'fieldName' from an exchange group named
!>  'groupName' within a domain's list of exchange groups.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_remove_field(domain, groupName, fieldName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      character (len=*), intent(in) :: fieldName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr, prevFieldListPtr
      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen, fieldNLen
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         nLen = len_trim(groupName)
         fieldNLen = len_trim(fieldName)
         DMPAR_DEBUG_WRITE(' -- Removing field ' // trim(fieldName) // ' from exchange group ' // trim(groupName))

         exchGroupPtr => domain % exchangeGroups
         do while (associated(exchGroupPtr))
            if ( nLen == exchGroupPtr % nLen) then
               if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
                  ! Remove after the transition to exchange_field_lists is complete
                  call mpas_pool_remove_config(exchGroupPtr % fieldPool, fieldName)

                  if ( associated(exchGroupPtr % fieldList) ) then
                     exchFieldListPtr => exchGroupPtr % fieldList
                     nullify(prevFieldListPtr)
                     do while ( associated(exchFieldListPtr) )
                        if ( fieldNLen == exchFieldListPtr % nLen ) then
                           if ( fieldName(1:fieldNLen) == exchFieldListPtr % fieldName(1:exchFieldListPtr % nLen) ) then
                              ! Remove exchFieldListPtr, as it's the field we're looking for.
                              if ( associated(prevFieldListPtr) ) then
                                 if ( associated(exchFieldListPtr % next) ) then
                                    prevFieldListPtr % next => exchFieldListPtr % next
                                 else
                                    nullify(prevFieldListPtr % next)
                                 end if
                              else
                                 if ( associated(exchFieldListPtr % next) ) then
                                    exchGroupPtr % fieldList => exchFieldListPtr % next
                                 else
                                    nullify(exchGroupPtr % fieldList)
                                 end if
                              end if

                              if ( associated(exchFieldListPtr % haloLayers) ) then
                                 deallocate(exchFieldListPtr % haloLayers)
                              end if
                              if ( associated(exchFieldListPtr % timeLevels) ) then
                                 deallocate(exchFieldListPtr % timeLevels)
                              end if
                              nullify(exchFieldListPtr % next)
                              deallocate(exchFieldListPtr)

                              ! Return, as we've successfully removed the field
                              ! from the field list.
                              return
                           end if
                        end if

                        prevFieldListPtr => exchFieldListPtr
                        exchFieldListPtr => exchFieldListPtr % next
                     end do
                  end if

                  return
               end if
            end if
            exchGroupPtr => exchGroupPtr % next
         end do

         write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
         write(stderrUnit, *) '       Cannot add field ''' // trim(fieldName) // ''' to group.'
         if ( present(iErr) ) then
            iErr = MPAS_DMPAR_MISSING_GROUP
         end if
      end if

   end subroutine mpas_dmpar_exch_group_remove_field!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_begin_halo_exch
!
!> \brief MPAS dmpar exchange group begin halo exchange routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine starts a halo exchange for an exchange group. This includes
!>  creating the buffers, packing the buffers, and starting the ISend / IRecv
!>  commands.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_begin_halo_exch(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      DMPAR_DEBUG_WRITE(' -- Trying to start halo exchange for group ' // trim(groupName))

      nLen = len_trim(groupName)

      exchGroupPtr => domain % exchangeGroups
      do while (associated(exchGroupPtr))
         if ( nLen == exchGroupPtr % nLen) then
            if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
               DMPAR_DEBUG_WRITE(' -- Starting halo exchange for group ' // trim(groupName))
               call mpas_dmpar_exch_group_build_buffers(domain % blocklist % allFields, exchGroupPtr)
               call mpas_dmpar_exch_group_start_recv(domain % dminfo, exchGroupPtr)
               call mpas_dmpar_exch_group_pack_buffers(domain % blocklist % allFields, exchGroupPtr)
               call mpas_dmpar_exch_group_start_send(domain % dminfo, exchGroupPtr)
               return
            end if
         end if
         exchGroupPtr => exchGroupPtr % next
      end do

      write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
      write(stderrUnit, *) '       Cannot perform halo exchange on group.'
      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_MISSING_GROUP
      end if

   end subroutine mpas_dmpar_exch_group_begin_halo_exch!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_halo_exch
!
!> \brief MPAS dmpar exchange group local halo exchange routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This routine performs the local exchange portion of a halo exchange using
!>  an exchange group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_halo_exch(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      nLen = len_trim(groupName)
      DMPAR_DEBUG_WRITE(' -- Trying to perform local copies for exchange group ' // trim(groupName))

      exchGroupPtr => domain % exchangeGroups
      do while (associated(exchGroupPtr))
         if ( nLen == exchGroupPtr % nLen) then
            if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
               DMPAR_DEBUG_WRITE(' -- Performing local copies for exchange group ' // trim(groupName))
               call mpas_dmpar_exch_group_local_exch_fields(domain % blocklist % allFields, exchGroupPtr)
               return
            end if
         end if
         exchGroupPtr => exchGroupPtr % next
      end do

      write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
      write(stderrUnit, *) '       Cannot perform halo exchange on group.'
      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_MISSING_GROUP
      end if

   end subroutine mpas_dmpar_exch_group_local_halo_exch!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_end_halo_exch
!
!> \brief MPAS dmpar exchange group end halo exchange routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This routine ends a halo exchange using an exchange group. This includes
!>  waiting for IRecv commands to complete for receive buffers, unpacking buffers,
!>  and deallocating buffers.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_end_halo_exch(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      nLen = len_trim(groupName)
      DMPAR_DEBUG_WRITE(' -- Trying to finish halo exchange for exchange group ' // trim(groupName))

      exchGroupPtr => domain % exchangeGroups
      do while (associated(exchGroupPtr))
         if ( nLen == exchGroupPtr % nLen) then
            if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
               DMPAR_DEBUG_WRITE(' -- Finishing halo exchange for exchange group ' // trim(groupName))
               call mpas_dmpar_exch_group_unpack_buffers(domain % blocklist % allFields, exchGroupPtr)
               return
            end if
         end if
         exchGroupPtr => exchGroupPtr % next
      end do

      write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
      write(stderrUnit, *) '       Cannot perform halo exchange on group.'
      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_MISSING_GROUP
      end if

   end subroutine mpas_dmpar_exch_group_end_halo_exch!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_full_halo_exch
!
!> \brief MPAS dmpar exchange group full halo exchange routine
!> \author Doug Jacobsen
!> \date   01/11/2016
!> \details
!>  This routine performs a full halo exchange on an exchange group.
!>  It is blocking, in that the routine doesn't return until the full
!>  exchange is complete.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_full_halo_exch(domain, groupName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: groupName
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      nLen = len_trim(groupName)
      DMPAR_DEBUG_WRITE(' -- Trying to perform a full exchange for group ' // trim(groupName))

      exchGroupPtr => domain % exchangeGroups
      do while (associated(exchGroupPtr))
         if ( nLen == exchGroupPtr % nLen) then
            if ( groupName(1:nLen) == exchGroupPtr % groupName(1:exchGroupPtr % nLen) ) then
               DMPAR_DEBUG_WRITE(' -- Performing a full exchange for group ' // trim(groupName))
               ! Setup exchange (should be the same as begin_halo_exch)
               call mpas_dmpar_exch_group_build_buffers(domain % blocklist % allFields, exchGroupPtr)
               call mpas_dmpar_exch_group_start_recv(domain % dminfo, exchGroupPtr)
               call mpas_dmpar_exch_group_pack_buffers(domain % blocklist % allFields, exchGroupPtr)
               call mpas_dmpar_exch_group_start_send(domain % dminfo, exchGroupPtr)

               ! Perform local copies (should be teh same as local_halo_exch)
               call mpas_dmpar_exch_group_local_exch_fields(domain % blocklist % allFields, exchGroupPtr)

               ! Finish the halo exchange (should be the same as end_halo_exch)
               call mpas_dmpar_exch_group_unpack_buffers(domain % blocklist % allFields, exchGroupPtr)

               ! Print out buffers for debugging
               !call mpas_dmpar_exch_group_print_buffers(exchGroupPtr)

               ! Destroy the buffers
               call mpas_dmpar_exch_group_destroy_buffers(exchGroupPtr)
               return
            end if
         end if
         exchGroupPtr => exchGroupPtr % next
      end do

      write(stderrUnit, *) 'ERROR: No exchange group found named ''' // trim(groupName) // '''.'
      write(stderrUnit, *) '       Cannot perform halo exchange on group.'
      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_MISSING_GROUP
      end if

   end subroutine mpas_dmpar_exch_group_full_halo_exch!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_field_halo_exch
!
!> \brief MPAS dmpar full halo exchange routine
!> \author Doug Jacobsen
!> \date   01/11/2016
!> \details
!>  This routine performs a full halo exchange on a specific field.
!>  It is blocking, in that the routine doesn't return until the full
!>  exchange is complete.
!>  It creates a temporary exchange group, adds this field to it, exchanges it,
!>  and destroys the group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_field_halo_exch(domain, fieldName, timeLevel, haloLayers, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: fieldName
      integer, optional, intent(in) :: timeLevel
      integer, dimension(:), optional, intent(in) :: haloLayers
      integer, optional, intent(out) :: iErr

      character (len=StrKIND) :: groupName
      type (mpas_exchange_group), pointer :: exchGroupPtr
      integer :: nLen

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      DMPAR_DEBUG_WRITE(' -- Performing a single field halo exchange for ' // trim(fieldName))

      groupName = 'TEMPSingleFieldGroup'
      call mpas_dmpar_exch_group_create(domain, groupName)

      call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName, timeLevel=timeLevel, haloLayers=haloLayers)

      call mpas_threading_barrier()

      call mpas_dmpar_exch_group_full_halo_exch(domain, groupName)

      call mpas_threading_barrier()

      call mpas_dmpar_exch_group_destroy(domain, groupName)

   end subroutine mpas_dmpar_field_halo_exch!}}}


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Private routines for performing the beginning of halo exchanges
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_build_buffers
!
!> \brief MPAS dmpar exchange group buffer construction routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine creates the buffers and communication lists for an exchange group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_build_buffers(allFieldsPool, exchangeGroup, iErr)!{{{
      type (mpas_pool_type), intent(in) :: allFieldsPool
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: intField1D
      type (field2DInteger), pointer :: intField2D
      type (field3DInteger), pointer :: intField3D
      type (field1DReal), pointer :: realField1D
      type (field2DReal), pointer :: realField2D
      type (field3DReal), pointer :: realField3D
      type (field4DReal), pointer :: realField4D
      type (field5DReal), pointer :: realField5D

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iHalo, iTimeLevel
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      ! Only master thread can create buffers.
      if ( threadNum == 0 ) then
         ! Allocate communication lists, and setup dead header nodes
         allocate(exchangeGroup % sendList)
         nullify(exchangeGroup % sendList % next)
         exchangeGroup % sendList % procID = -1
         exchangeGroup % sendList % nList = 0

         allocate(exchangeGroup % recvList)
         nullify(exchangeGroup % recvList % next)
         exchangeGroup % recvList % procID = -1
         exchangeGroup % recvList % nList = 0

         exchFieldListPtr => exchangeGroup % fieldList
         do while ( associated(exchFieldListPtr) )
            DMPAR_DEBUG_WRITE(' -- Building buffers for field ' // trim(exchFieldListPtr % fieldName))
            do iTimeLevel = 1, size(exchFieldListPtr % timeLevels)
               if ( exchFieldListPtr % timeLevels(iTimeLevel) ) then
                  if ( exchFieldListPtr % fieldType == MPAS_POOL_REAL ) then
                     if ( exchFieldListPtr % nDims == 1 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField1D, iTimeLevel)

                        if ( realField1D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField1D % sendList, exchangeGroup % sendList, &
                                                                                realField1D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField1D % recvList, exchangeGroup % recvList, &
                                                                                realField1D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 2 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField2D, iTimeLevel)

                        if ( realField2D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField2D % sendList, exchangeGroup % sendList, &
                                                                                realField2D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField2D % recvList, exchangeGroup % recvList, &
                                                                                realField2D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 3 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField3D, iTimeLevel)

                        if ( realField3D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField3D % sendList, exchangeGroup % sendList, &
                                                                                realField3D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField3D % recvList, exchangeGroup % recvList, &
                                                                                realField3D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 4 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField4D, iTimeLevel)

                        if ( realField4D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField4D % sendList, exchangeGroup % sendList, &
                                                                                realField4D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField4D % recvList, exchangeGroup % recvList, &
                                                                                realField4D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 5 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField5D, iTimeLevel)

                        if ( realField5D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField5D % sendList, exchangeGroup % sendList, &
                                                                                realField5D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(realField5D % recvList, exchangeGroup % recvList, &
                                                                                realField5D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     end if
                  else if ( exchFieldListPtr % fieldType == MPAS_POOL_INTEGER ) then
                     if ( exchFieldListPtr % nDims == 1 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField1D, iTimeLevel)

                        if ( intField1D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField1D % sendList, exchangeGroup % sendList, &
                                                                                intField1D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField1D % recvList, exchangeGroup % recvList, &
                                                                                intField1D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 2 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField2D, iTimeLevel)

                        if ( intField2D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField2D % sendList, exchangeGroup % sendList, &
                                                                                intField2D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField2D % recvList, exchangeGroup % recvList, &
                                                                                intField2D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     else if ( exchFieldListPtr % nDims == 3 ) then
                        call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField3D, iTimeLevel)

                        if ( intField3D % isActive ) then
                           do iHalo = 1, size(exchFieldListPtr % haloLayers)
                              if ( exchFieldListPtr % haloLayers(iHalo) ) then
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField3D % sendList, exchangeGroup % sendList, &
                                                                                intField3D % dimSizes, iHalo)
                                 call mpas_dmpar_exch_group_aggregate_exch_list(intField3D % recvList, exchangeGroup % recvList, &
                                                                                intField3D % dimSizes, iHalo, .true.)
                              end if
                           end do
                        end if
                     end if
                  end if
               end if
            end do

            exchFieldListPtr => exchFieldListPtr % next
         end do

         ! Remove dead header nodes on communication lists
         commListPtr => exchangeGroup % sendList
         if ( associated(exchangeGroup % sendList % next) ) then
            exchangeGroup % sendList => exchangeGroup % sendList % next
         else
            nullify(exchangeGroup % sendList)
         end if
         deallocate(commListPtr)

         commListPtr => exchangeGroup % recvList
         if ( associated(exchangeGroup % recvList % next) ) then
            exchangeGroup % recvList => exchangeGroup % recvList % next
         else
            nullify(exchangeGroup % recvList)
         end if
         deallocate(commListPtr)


         ! Allocate buffers for each processor's communication list
         commListPtr => exchangeGroup % sendList
         do while ( associated(commListPtr) )
            if ( associated(commListPtr % rbuffer) ) then
               write(stderrUnit, *) 'ERROR: Communication buffer already exists. A halo exchange might be ' // &
                                    'in progress for group ''' // trim(exchangeGroup % groupName) // '''.'
               if ( present(iErr) ) then
                  iErr = MPAS_DMPAR_BUFFER_EXISTS
               end if
               return
            end if
            allocate(commListPtr % rbuffer(commListPtr % nList))
            commListPtr % bufferOffset = 0

            commListPtr => commListPtr % next
         end do

         commListPtr => exchangeGroup % recvList
         do while ( associated(commListPtr) )
            if ( associated(commListPtr % rbuffer) ) then
               write(stderrUnit, *) 'ERROR: Communication buffer already exists. A halo exchange might be ' // &
                                    'in progress for group ''' // trim(exchangeGroup % groupName) // '''.'
               if ( present(iErr) ) then
                  iErr = MPAS_DMPAR_BUFFER_EXISTS
               end if
               return
            end if
            allocate(commListPtr % rbuffer(commListPtr % nList))
            commListPtr % bufferOffset = 0

            commListPtr => commListPtr % next
         end do
      end if

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_build_buffers!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_aggregate_exch_list
!
!> \brief MPAS dmpar exchange group exchange list aggregation routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine aggregates a multihalo exchange list into a communication list.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_aggregate_exch_list(multiHaloExchList, commList, dimSizes, haloLayer, recvList, iErr)!{{{
      type (mpas_multihalo_exchange_list), pointer :: multiHaloExchList
      type (mpas_communication_list), pointer :: commList
      integer,  dimension(:), intent(in) :: dimSizes
      integer, intent(in) :: haloLayer
      logical, optional, intent(in) :: recvList
      integer, optional, intent(out) :: iErr

      type (mpas_multihalo_exchange_list), pointer :: multiHaloExchListPtr
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr, commListPtr2

      logical :: comm_list_found = .false.
      logical :: recvListLocal
      integer :: bufferOffset, nAdded
      integer :: dimProd

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      if ( present(recvList) ) then
         recvListLocal = recvList
      else
         recvListLocal = .false.
      end if

      call mpas_dmpar_build_dim_size_product(dimSizes, dimProd)

      multiHaloExchListPtr => multiHaloExchList
      do while ( associated(multiHaloExchListPtr) )
         ! Aggregate this fields send list for this halo into the buffer size
         exchListPtr => multiHaloExchListPtr % halos(haloLayer) % exchList
         do while (associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => commList
            do while(associated(commListPtr))
               if ( commListPtr % procID == exchListPtr % endPointID ) then
                  comm_list_found = .true.
                  if ( .not. recvListLocal ) then
                     commListPtr % nList = commListPtr % nList + exchListPtr % nList * dimProd
                  end if
                  exit
               end if
               commListPtr => commListPtr % next
            end do

            ! Set up a new communication list for this exchange list
            ! because it was not found. Put it at the head of the list to make insertion quick.
            if ( .not. comm_list_found ) then
               commListPtr => commList

               allocate(commListPtr2)
               if ( associated(commListPtr % next) ) then
                  commListPtr2 % next => commListPtr % next
               else
                  nullify(commListPtr2 % next)
               end if
               commListPtr2 % procID = exchListPtr % endPointID
               if ( .not. recvListLocal ) then
                  commListPtr2 % nList = exchListPtr % nList * dimProd
               else
                  commListPtr2 % nList = 0
               end if
               nullify(commListPtr2 % rbuffer)
               nullify(commListPtr2 % ibuffer)

               commListPtr % next => commListPtr2
            end if

            exchListPtr => exchListPtr % next
         end do
         multiHaloExchListPtr => multiHaloExchListPtr % next
      end do

      ! If this is a recive list, nlist needs to be build differently.
      ! Determine size of receive lists
      if ( recvListLocal ) then
         commListPtr => commList
         do while(associated(commListPtr))
            bufferOffset = 0
            nAdded = 0

            multiHaloExchListPtr => multiHaloExchList
            do while(associated(multiHaloExchListPtr))
               exchListPtr => multiHaloExchListPtr % halos(haloLayer) % exchList
               do while(associated(exchListPtr))
                  if(exchListPtr % endPointID == commListPtr % procID) then
                     nAdded = max(nAdded, maxval(exchListPtr % srcList) * dimProd)
                  end if
                  exchListPtr => exchListPtr % next
               end do

               multiHaloExchListPtr => multiHaloExchListPtr % next
            end do
            bufferOffset = bufferOffset + nAdded
            commListPtr % nList = commListPtr % nList + nAdded
            commListPtr % bufferOffset = 0

            commListPtr => commListPtr % next
         end do  ! commListPtr
      end if


   end subroutine mpas_dmpar_exch_group_aggregate_exch_list!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_build_dim_size_product
!
!> \brief MPAS dmpar build dimension size product routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine determines the dimension size product for a field,
!>  which is used when building the size of a halo exchange buffer.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_build_dim_size_product(dimSizes, dimSizeProd)!{{{
      integer, dimension(:), intent(in) :: dimSizes
      integer, intent(out) :: dimSizeProd

      integer :: iDimen

      dimSizeProd = 1

      do iDimen = 1, size(dimSizes) - 1
         dimSizeProd = dimSizeProd * dimSizes(iDimen)
      end do

   end subroutine mpas_dmpar_build_dim_size_product!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_start_recv
!
!> \brief MPAS dmpar exchange group start irecv routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine starts the irecv commands for each communication in an exchange group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_start_recv(dminfo, exchangeGroup)!{{{
      type (dm_info), intent(in) :: dminfo
      type (mpas_exchange_group), intent(inout) :: exchangeGroup

      type (mpas_communication_list), pointer :: commListPtr

      integer :: threadNum
      integer :: mpi_ierr

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         commListPtr => exchangeGroup % recvList
         do while ( associated(commListPtr) )
            DMPAR_DEBUG_WRITE('    -- Starting recv: ' COMMA commListPtr % procID COMMA commListPtr % nList COMMA size(commListPtr % rbuffer) )
            call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, &
                           dminfo % comm, commListPtr % reqID, mpi_ierr)

            commListPtr => commListPtr % next
         end do
      end if

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_start_recv!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_start_send
!
!> \brief MPAS dmpar exchange group start isend routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine starts the isend commands for each communication in an exchange group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_start_send(dminfo, exchangeGroup)!{{{
      type (dm_info), intent(in) :: dminfo
      type (mpas_exchange_group), intent(inout) :: exchangeGroup

      type (mpas_communication_list), pointer :: commListPtr

      integer :: threadNum
      integer :: mpi_ierr

      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         commListPtr => exchangeGroup % sendList
         do while ( associated(commListPtr) )
            DMPAR_DEBUG_WRITE('    -- Starting send: ' COMMA commListPtr % procID COMMA commListPtr % nList COMMA size(commListPtr % rbuffer) )
            call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, &
                           dminfo % comm, commListPtr % reqID, mpi_ierr)

            commListPtr => commListPtr % next
         end do
      end if

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_start_send!}}}


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Private routines for group operations on a list of fields
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffers
!
!> \brief MPAS dmpar exchange group buffer pack routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs all fields into their send buffers for an exchange group
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffers(allFieldsPool, exchangeGroup, iErr)!{{{
      type (mpas_pool_type), intent(in) :: allFieldsPool
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      integer, optional, intent(out) :: iErr

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr

      type (field1DInteger), pointer :: intField1D
      type (field2DInteger), pointer :: intField2D
      type (field3DInteger), pointer :: intField3D
      type (field1DReal), pointer :: realField1D
      type (field2DReal), pointer :: realField2D
      type (field3DReal), pointer :: realField3D
      type (field4DReal), pointer :: realField4D
      type (field5DReal), pointer :: realField5D

      integer :: iTimeLevel, iHalo

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      ! Pack each feild entirely before moving on to the next field
      exchFieldListPtr => exchangeGroup % fieldList
      do while ( associated(exchFieldListPtr) )
         DMPAR_DEBUG_WRITE(' -- Packing buffers for field ' // trim(exchFieldListPtr % fieldName))
         do iTimeLevel = 1, size(exchFieldListPtr % timeLevels)
            if ( exchFieldListPtr % timeLevels(iTimeLevel) ) then

               if ( exchFieldListPtr % fieldType == MPAS_POOL_REAL ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField1D, iTimeLevel)
                     if ( realField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, realField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField2D, iTimeLevel)
                     if ( realField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, realField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField3D, iTimeLevel)
                     if ( realField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, realField3D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 4 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField4D, iTimeLevel)
                     if ( realField4D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, realField4D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 5 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField5D, iTimeLevel)
                     if ( realField5D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, realField5D, iHalo)
                           end if
                        end do
                     end if
                  end if
               else if ( exchFieldListPtr % fieldType == MPAS_POOL_INTEGER ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField1D, iTimeLevel)
                     if ( intField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, intField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField2D, iTimeLevel)
                     if ( intField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, intField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField3D, iTimeLevel)
                     if ( intField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_pack_buffer_field(exchangeGroup, intField3D, iHalo)
                           end if
                        end do
                     end if
                  end if
               end if
            end if
         end do

         exchFieldListPtr => exchFieldListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffers!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_fields
!
!> \brief MPAS dmpar exchange group local exchange routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine performs the actual local exchanges for each field in
!>  an exchange group.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_fields(allFieldsPool, exchangeGroup, iErr)!{{{
      type (mpas_pool_type), intent(in) :: allFieldsPool
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: intField1D
      type (field2DInteger), pointer :: intField2D
      type (field3DInteger), pointer :: intField3D
      type (field1DReal), pointer :: realField1D
      type (field2DReal), pointer :: realField2D
      type (field3DReal), pointer :: realField3D
      type (field4DReal), pointer :: realField4D
      type (field5DReal), pointer :: realField5D

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr

      integer :: iTimeLevel, iHalo

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      exchFieldListPtr => exchangeGroup % fieldList
      do while ( associated(exchFieldListPtr) )
         DMPAR_DEBUG_WRITE(' -- Local copy for field ' // trim(exchFieldListPtr % fieldName))
         do iTimeLevel = 1, size(exchFieldListPtr % timeLevels)
            if ( exchFieldListPtr % timeLevels(iTimeLevel) ) then
               if ( exchFieldListPtr % fieldType == MPAS_POOL_REAL ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField1D, iTimeLevel)
                     if ( realField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(realField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField2D, iTimeLevel)
                     if ( realField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(realField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField3D, iTimeLevel)
                     if ( realField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(realField3D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 4 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField4D, iTimeLevel)
                     if ( realField4D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(realField4D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 5 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField5D, iTimeLevel)
                     if ( realField5D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(realField5D, iHalo)
                           end if
                        end do
                     end if
                  end if
               else if ( exchFieldListPtr % fieldType == MPAS_POOL_INTEGER ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField1D, iTimeLevel)
                     if ( intField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(intField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField2D, iTimeLevel)
                     if ( intField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(intField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField3D, iTimeLevel)
                     if ( intField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_local_exch_field(intField3D, iHalo)
                           end if
                        end do
                     end if
                  end if
               end if
            end if
         end do

         exchFieldListPtr => exchFieldListPtr % next
      end do

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_local_exch_fields!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffers
!
!> \brief MPAS dmpar exchange group buffer unpack routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks buffers. Additionally, it includes MPI_Wait
!>  commands to finish receiving messages before unpacking buffers.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffers(allFieldsPool, exchangeGroup, iErr)!{{{
      type (mpas_pool_type), intent(in) :: allFieldsPool
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: intField1D
      type (field2DInteger), pointer :: intField2D
      type (field3DInteger), pointer :: intField3D
      type (field1DReal), pointer :: realField1D
      type (field2DReal), pointer :: realField2D
      type (field3DReal), pointer :: realField3D
      type (field4DReal), pointer :: realField4D
      type (field5DReal), pointer :: realField5D

      type (mpas_exchange_field_list), pointer :: exchFieldListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: mpi_ierr
      integer :: threadNum
      integer :: iTimeLevel, iHalo

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      exchFieldListPtr => exchangeGroup % fieldList
      do while ( associated(exchFieldListPtr) )
         DMPAR_DEBUG_WRITE(' -- Unpacking buffers for field ' // trim(exchFieldListPtr % fieldName))
         do iTimeLevel = 1, size(exchFieldListPtr % timeLevels)
            if ( exchFieldListPtr % timeLevels(iTimeLevel) ) then
               if ( exchFieldListPtr % fieldType == MPAS_POOL_REAL ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField1D, iTimeLevel)
                     if ( realField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, realField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField2D, iTimeLevel)
                     if ( realField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, realField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField3D, iTimeLevel)
                     if ( realField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, realField3D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 4 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField4D, iTimeLevel)
                     if ( realField4D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, realField4D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 5 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, realField5D, iTimeLevel)
                     if ( realField5D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, realField5D, iHalo)
                           end if
                        end do
                     end if
                  end if
               else if ( exchFieldListPtr % fieldType == MPAS_POOL_INTEGER ) then
                  if ( exchFieldListPtr % nDims == 1 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField1D, iTimeLevel)
                     if ( intField1D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, intField1D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 2 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField2D, iTimeLevel)
                     if ( intField2D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, intField2D, iHalo)
                           end if
                        end do
                     end if
                  else if ( exchFieldListPtr % nDims == 3 ) then
                     call mpas_pool_get_field(allFieldsPool, exchFieldListPtr % fieldName, intField3D, iTimeLevel)
                     if ( intField3D % isActive ) then
                        do iHalo = 1, size(exchFieldListPtr % haloLayers)
                           if ( exchFieldListPtr % haloLayers(iHalo) ) then
                              call mpas_dmpar_exch_group_unpack_buffer_field(exchangeGroup, intField3D, iHalo)
                           end if
                        end do
                     end if
                  end if
               end if
            end if
         end do


         exchFieldListPtr => exchFieldListPtr % next
      end do

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_unpack_buffers!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_destroy_buffers
!
!> \brief MPAS dmpar exchange group destroy buffers routine
!> \author Doug Jacobsen
!> \date   04/07/2016
!> \details
!>  This private routine destroys buffers. Additionally, it includes MPI_Wait
!>  commands to finish receiving messages before destroying buffers.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_destroy_buffers(exchangeGroup, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      integer, optional, intent(out) :: iErr

      type (mpas_communication_list), pointer :: commListPtr

      integer :: mpi_ierr
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      call mpas_threading_barrier()

      if ( threadNum == 0 ) then
         ! Wait for isends to finish
         commListPtr => exchangeGroup % sendList
         do while ( associated(commListPtr) )
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

            commListPtr => commListPtr % next
         end do

         ! Destroy communication lists
         call mpas_dmpar_destroy_communication_list(exchangeGroup % sendList)
         call mpas_dmpar_destroy_communication_list(exchangeGroup % recvList)
      end if

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_destroy_buffers!}}}



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Private routines for packing send buffers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field1d_integer
!
!> \brief MPAS dmpar exchange pack 1D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 1D integer field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field1d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field1DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(iBuffer)
                  do iExch = 1, exchListPtr % nList
                     iBuffer = exchListPtr % destList(iExch) + bufferOffset
                     commListPtr % rbuffer(iBuffer) = transfer(fieldCursor % array(exchListPtr % srcList(iExch)), &
                                                               commListPtr % rbuffer(1))
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field1d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field2d_integer
!
!> \brief MPAS dmpar exchange pack 2D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 2D integer field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field2d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field2DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(1)
                        iBuffer = (exchListPtr % destList(iExch)-1) * fieldCursor % dimSizes(1) + j + bufferOffset
                        commListPtr % rbuffer(iBuffer) = transfer(fieldCursor % array(j, exchListPtr % srcList(iExch)), &
                                                                  commListPtr % rbuffer(1))
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field2d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field3d_integer
!
!> \brief MPAS dmpar exchange pack 3D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 3D integer field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field3d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field3DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j, k

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(2)
                        do k = 1, fieldCursor % dimSizes(1)
                           iBuffer = (exchListPtr % destList(iExch) - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                   + (j - 1) * fieldCursor % dimSizes(1) + k + bufferOffset
                           commListPtr % rbuffer(iBuffer) = transfer(fieldCursor % array(k, j, exchListPtr % srcList(iExch)), &
                                                                     commListPtr % rbuffer(1))
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field3d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field1d_real
!
!> \brief MPAS dmpar exchange pack 1D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 1D real field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field1d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field1DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(iBuffer)
                  do iExch = 1, exchListPtr % nList
                     iBuffer = exchListPtr % destList(iExch) + bufferOffset
                     commListPtr % rbuffer(iBuffer) = fieldCursor % array(exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field1d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field2d_real
!
!> \brief MPAS dmpar exchange pack 2D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 2D real field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field2d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field2DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(1)
                        iBuffer = (exchListPtr % destList(iExch)-1) * fieldCursor % dimSizes(1) + j + bufferOffset
                        commListPtr % rbuffer(iBuffer) = fieldCursor % array(j, exchListPtr % srcList(iExch))
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field2d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field3d_real
!
!> \brief MPAS dmpar exchange pack 3D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 3D real field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field3d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field3DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j, k

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(2)
                        do k = 1, fieldCursor % dimSizes(1)
                           iBuffer = (exchListPtr % destList(iExch) - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                   + (j - 1) * fieldCursor % dimSizes(1) + k + bufferOffset
                           commListPtr % rbuffer(iBuffer) = fieldCursor % array(k, j, exchListPtr % srcList(iExch))
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field3d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field4d_real
!
!> \brief MPAS dmpar exchange pack 4D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 4D real field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field4d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field4DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field4DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j, k, l

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, l, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(3)
                        do k = 1, fieldCursor % dimSizes(2)
                           do l = 1, fieldCursor % dimSizes(1)
                              iBuffer = (exchListPtr % destList(iExch) - 1) * fieldCursor % dimSizes(1) &
                                      * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                      + (j - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                      + (k - 1) * fieldCursor % dimSizes(1) + l + bufferOffset
                              commListPtr % rbuffer(iBuffer) = fieldCursor % array(l, k, j, exchListPtr % srcList(iExch))
                           end do
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                         * fieldCursor % dimSizes(3)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field4d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_pack_buffer_field5d_real
!
!> \brief MPAS dmpar exchange pack 5D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine packs a 5D real field into the send buffers
!>  It packs a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_pack_buffer_field5d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field5DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field5DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded
      integer :: j, k, l, m

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      commListPtr => exchangeGroup % sendList
      do while ( associated(commListPtr) )
         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, l, m, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(4)
                        do k = 1, fieldCursor % dimSizes(3)
                           do l = 1, fieldCursor % dimSizes(2)
                              do m = 1, fieldCursor % dimSizes(1)
                                 iBuffer = (exchListPtr % destList(iExch) - 1) * fieldCursor % dimSizes(1) &
                                         * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4)&
                                         + (j - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                   * fieldCursor % dimSizes(3) &
                                         + (k - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                         + (l - 1) * fieldCursor % dimSizes(1) + m + bufferOffset
                                 commListPtr % rbuffer(iBuffer) = fieldCursor % array(m, l, k, j, exchListPtr % srcList(iExch))
                              end do
                           end do
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = nAdded + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                         * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4)
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         !$omp master
         commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         !$omp end master
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_pack_buffer_field5d_real!}}}


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Private routines for performing local exchanges
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field1d_integer
!
!> \brief MPAS dmpar exchange pack 1D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 1D integer field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field1d_integer(field, haloLayer, iErr)!{{{
      type (field1DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(exchListPtr % destList(iExch)) = fieldCursor % array(exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field1d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field2d_integer
!
!> \brief MPAS dmpar exchange pack 2D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 2D integer field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field2d_integer(field, haloLayer, iErr)!{{{
      type (field2DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DInteger), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, exchListPtr % destList(iExch)) = fieldCursor % array(:, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field2d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field3d_integer
!
!> \brief MPAS dmpar exchange pack 3D integer field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 3D integer field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field3d_integer(field, haloLayer, iErr)!{{{
      type (field3DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DInteger), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, :, exchListPtr % destList(iExch)) = &
                                  fieldCursor % array(:, :, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field3d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field1d_real
!
!> \brief MPAS dmpar exchange pack 1D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 1D real field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field1d_real(field, haloLayer, iErr)!{{{
      type (field1DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DReal), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(exchListPtr % destList(iExch)) = fieldCursor % array(exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field1d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field2d_real
!
!> \brief MPAS dmpar exchange pack 2D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 2D real field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field2d_real(field, haloLayer, iErr)!{{{
      type (field2DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DReal), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, exchListPtr % destList(iExch)) = fieldCursor % array(:, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field2d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field3d_real
!
!> \brief MPAS dmpar exchange pack 3D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 3D real field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field3d_real(field, haloLayer, iErr)!{{{
      type (field3DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DReal), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, :, exchListPtr % destList(iExch)) = &
                                  fieldCursor % array(:, :, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field3d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field4d_real
!
!> \brief MPAS dmpar exchange pack 4D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 4D real field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field4d_real(field, haloLayer, iErr)!{{{
      type (field4DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field4DReal), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, :, :, exchListPtr % destList(iExch)) = &
                                  fieldCursor % array(:, :, :, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field4d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_local_exch_field5d_real
!
!> \brief MPAS dmpar exchange pack 5D real field routine
!> \author Doug Jacobsen
!> \date   01/05/2016
!> \details
!>  This private routine exchanges a 5D real field between local blocks
!>  It exchanges a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_local_exch_field5d_real(field, haloLayer, iErr)!{{{
      type (field5DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field5DReal), pointer :: fieldCursor, fieldCursor2

      type (mpas_exchange_list), pointer :: exchListPtr

      integer :: iExch

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      fieldCursor => field
      do while ( associated(fieldCursor) )
         exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList

         do while ( associated(exchListPtr) )
            fieldCursor2 => field
            do while ( associated(fieldCursor2) )
               if ( exchListPtr % endPointID == fieldCursor2 % block % localBlockID ) then
                  !$omp do schedule(runtime)
                  do iExch = 1, exchListPtr % nList
                     fieldCursor2 % array(:, :, :, :, exchListPtr % destList(iExch)) = &
                                  fieldCursor % array(:, :, :, :, exchListPtr % srcList(iExch))
                  end do
                  !$omp end do
               end if

               fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
         end do

         call mpas_threading_barrier()

         fieldCursor => fieldCursor % next
      end do

   end subroutine mpas_dmpar_exch_group_local_exch_field5d_real!}}}


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Private routines for unpacking receive buffers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field1d_integer
!
!> \brief MPAS dmpar exchange unpack 1D integer field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 1D integer field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field1d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field1DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(iBuffer)
                  do iExch = 1, exchListPtr % nList
                     iBuffer = exchListPtr % srcList(iExch) + bufferOffset
                     fieldCursor % array(exchListPtr % destList(iExch)) = transfer(commListPtr % rbuffer(iBuffer), &
                                                                                   fieldCursor % array(1))
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field1d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field2d_integer
!
!> \brief MPAS dmpar exchange unpack 2D integer field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 2D integer field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field2d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field2DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(1)
                        iBuffer = (exchListPtr % srcList(iExch)-1) * fieldCursor % dimSizes(1) + j + bufferOffset
                        fieldCursor % array(j, exchListPtr % destList(iExch)) = transfer(commListPtr % rbuffer(iBuffer), &
                                                                                         fieldCursor % array(1,1))
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field2d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field3d_integer
!
!> \brief MPAS dmpar exchange unpack 3D integer field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 3D integer field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field3d_integer(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field3DInteger), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DInteger), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j, k
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(2)
                        do k = 1, fieldCursor % dimSizes(1)
                           iBuffer = (exchListPtr % srcList(iExch) - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                   + (j - 1) * fieldCursor % dimSizes(1) + k + bufferOffset
                           fieldCursor % array(k, j, exchListPtr % destList(iExch)) = transfer(commListPtr % rbuffer(iBuffer), &
                                                                                               fieldCursor % array(1,1,1))
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field3d_integer!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field1d_real
!
!> \brief MPAS dmpar exchange unpack 1D real field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 1D real field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field1d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field1DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field1DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(iBuffer)
                  do iExch = 1, exchListPtr % nList
                     iBuffer = exchListPtr % srcList(iExch) + bufferOffset
                     fieldCursor % array(exchListPtr % destList(iExch)) = commListPtr % rbuffer(iBuffer)
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field1d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field2d_real
!
!> \brief MPAS dmpar exchange unpack 2D real field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 2D real field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field2d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field2DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field2DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(1)
                        iBuffer = (exchListPtr % srcList(iExch)-1) * fieldCursor % dimSizes(1) + j + bufferOffset
                        fieldCursor % array(j, exchListPtr % destList(iExch)) = commListPtr % rbuffer(iBuffer)
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field2d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field3d_real
!
!> \brief MPAS dmpar exchange unpack 3D real field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 3D real field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field3d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field3DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field3DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j, k
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(2)
                        do k = 1, fieldCursor % dimSizes(1)
                           iBuffer = (exchListPtr % srcList(iExch) - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                   + (j - 1) * fieldCursor % dimSizes(1) + k + bufferOffset
                           fieldCursor % array(k, j, exchListPtr % destList(iExch)) = commListPtr % rbuffer(iBuffer)
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field3d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field4d_real
!
!> \brief MPAS dmpar exchange unpack 4D real field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 4D real field from the receive buffers
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field4d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field4DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field4DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j, k, l
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, l, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(3)
                        do k = 1, fieldCursor % dimSizes(2)
                           do l = 1, fieldCursor % dimSizes(1)
                              iBuffer = (exchListPtr % srcList(iExch) - 1) * fieldCursor % dimSizes(1) &
                                      * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                      + (j - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                      + (k - 1) * fieldCursor % dimSizes(1) + l + bufferOffset
                              fieldCursor % array(l, k, j, exchListPtr % destList(iExch)) = commListPtr % rbuffer(iBuffer)
                           end do
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) &
                         * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field4d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_unpack_buffer_field5d_real
!
!> \brief MPAS dmpar exchange unpack 5D real field routine
!> \author Doug Jacobsen
!> \date   01/06/2016
!> \details
!>  This private routine unpacks a 5D real field from the receive buffers.
!>  It unpacks a specific halo layer, as defined by the haloLayer input argument.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_unpack_buffer_field5d_real(exchangeGroup, field, haloLayer, iErr)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup
      type (field5DReal), pointer :: field
      integer, intent(in) :: haloLayer
      integer, optional, intent(out) :: iErr

      type (field5DReal), pointer :: fieldCursor

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iExch, iBuffer, bufferOffset
      integer :: nAdded, mpi_ierr
      integer :: j, k, l, m
      integer :: threadNum

      if ( present(iErr) ) then
         iErr = MPAS_DMPAR_NOERR
      end if

      threadNum = mpas_threading_get_thread_num()

      commListPtr => exchangeGroup % recvList
      do while ( associated(commListPtr) )
         if ( threadNum == 0 ) then
            call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
         end if
         call mpas_threading_barrier()

         bufferOffset = commListPtr % bufferOffset
         nAdded = 0
         fieldCursor => field
         do while ( associated(fieldCursor) )
            exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
            do while ( associated(exchListPtr) )
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  !$omp do schedule(runtime) private(j, k, l, m, iBuffer)
                  do iExch = 1, exchListPtr % nList
                     do j = 1, fieldCursor % dimSizes(4)
                        do k = 1, fieldCursor % dimSizes(3)
                           do l = 1, fieldCursor % dimSizes(2)
                              do m = 1, fieldCursor % dimSizes(1)
                                 iBuffer = (exchListPtr % srcList(iExch) - 1) * fieldCursor % dimSizes(1) &
                                         * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4)&
                                         + (j - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                   * fieldCursor % dimSizes(3) &
                                         + (k - 1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                         + (l - 1) * fieldCursor % dimSizes(1) + m + bufferOffset
                                 fieldCursor % array(m, l, k, j, exchListPtr % destList(iExch)) = commListPtr % rbuffer(iBuffer)
                              end do
                           end do
                        end do
                     end do
                  end do
                  !$omp end do
                  nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                         * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4))
               end if
               exchListPtr => exchListPtr % next
            end do
            fieldCursor => fieldCursor % next
         end do
         bufferOffset = bufferOffset + nAdded
         if ( threadNum == 0 ) then
            commListPtr % bufferOffset = commListPtr % bufferOffset + nAdded
         end if
         call mpas_threading_barrier()
         commListPtr => commListPtr % next
      end do

   end subroutine mpas_dmpar_exch_group_unpack_buffer_field5d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_dmpar_exch_group_print_buffers
!
!> \brief MPAS dmpar exchange group buffer info routine.
!> \author Doug Jacobsen
!> \date   04/07/2016
!> \details
!>  This private routine prints out information about communication list buffers.
!>  It is mostly used for debugging.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_exch_group_print_buffers(exchangeGroup)!{{{
      type (mpas_exchange_group), intent(inout) :: exchangeGroup

      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: commListPtr

      integer :: iBuffer, threadNum

      threadNum = mpas_threading_get_thread_num()

      call mpas_threading_barrier()

      if ( threadNum == 0 ) then
         write(stderrUnit, *) ''
         write(stderrUnit, *) ' -------- Buffer Summary -------'
         write(stderrUnit, *) ' Thread: ', threadNum, ' writing out send lists'
         commListPtr => exchangeGroup % sendList
         do while ( associated(commListPtr) )
            write(stderrUnit, *) ''
            write(stderrUnit, *) '    Send list header: '
            write(stderrUnit, *) '         proc: ', commListPtr % procID
            write(stderrUnit, *) '         size check: ', commListPtr % nlist, size( commListPtr % rbuffer )
            write(stderrUnit, *) '         bufferOffset: ', commListPtr % bufferOffset
            write(stderrUnit, *) '         reqId: ', commListPtr % reqId
            write(stderrUnit, *) '         ibuffer assc: ', associated( commListPtr % ibuffer )
            write(stderrUnit, *) '         rbuffer assc: ', associated( commListPtr % rbuffer )
            write(stderrUnit, *) '         next assc: ', associated( commListPtr % next )
            do iBuffer = 1, size(commListPtr % rbuffer)
               write(stderrUnit, *) '   IDX: ', iBuffer, ' VAL: ', commListPtr % rbuffer(iBuffer)
            end do
            commListPtr => commListPtr % next
         end do

         write(stderrUnit, *) ''
         write(stderrUnit, *) ' Thread: ', threadNum, ' writing out recv lists'
         commListPtr => exchangeGroup % recvList
         do while ( associated(commListPtr) )
            write(stderrUnit, *) ''
            write(stderrUnit, *) '    Recv list header: '
            write(stderrUnit, *) '         proc: ', commListPtr % procID
            write(stderrUnit, *) '         size check: ', commListPtr % nlist, size( commListPtr % rbuffer )
            write(stderrUnit, *) '         bufferOffset: ', commListPtr % bufferOffset
            write(stderrUnit, *) '         reqId: ', commListPtr % reqId
            write(stderrUnit, *) '         ibuffer assc: ', associated( commListPtr % ibuffer )
            write(stderrUnit, *) '         rbuffer assc: ', associated( commListPtr % rbuffer )
            write(stderrUnit, *) '         next assc: ', associated( commListPtr % next )
            do iBuffer = 1, size(commListPtr % rbuffer)
               write(stderrUnit, *) '   IDX: ', iBuffer, ' VAL: ', commListPtr % rbuffer(iBuffer)
            end do
            commListPtr => commListPtr % next
         end do

         write(stderrUnit, *) ''
         write(stderrUnit, *) ' -------- End Buffer Summary -------'
         write(stderrUnit, *) ''

      end if

      call mpas_threading_barrier()

   end subroutine mpas_dmpar_exch_group_print_buffers!}}}


end module mpas_dmpar


!-----------------------------------------------------------------------
!  routine mpas_dmpar_global_abort
!
!> \brief MPAS dmpar global abort routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine aborts MPI. A call to it kills the model through the use of
!>  MPI_Abort on the world communicator, and outputs a message.
!
!-----------------------------------------------------------------------
subroutine mpas_dmpar_global_abort(mesg)!{{{

   use mpas_kind_types, only : StrKIND
   use mpas_io_units, only : stderrUnit, stdoutUnit
   use mpas_threading, only : mpas_threading_get_thread_num

#ifdef _MPI
#ifndef NOMPIMOD
   use mpi
#endif
#endif

   implicit none

#ifdef _MPI
#ifdef NOMPIMOD
include 'mpif.h'
#endif
#endif

   character (len=*), intent(in) :: mesg !< Input: Abort message
   integer :: threadNum

#ifdef _MPI
   integer :: mpi_ierr, mpi_errcode, my_proc_id
#endif

   character(len=StrKIND) :: errorfile

   threadNum = mpas_threading_get_thread_num()

   if ( threadNum == 0 ) then
     write(stderrUnit,*) trim(mesg)
   end if

   ! close log units
   close(stderrUnit)
   close(stdoutUnit)

   ! write out error message to critical error message file
#ifdef _MPI
   call MPI_COMM_RANK(MPI_COMM_WORLD, my_proc_id, mpi_ierr)
   write(errorfile,fmt='(a,i6.6,a)') 'mpas_critical_error.', my_proc_id, '.log'
#else
   errorfile = 'mpas_critical_error.log'
#endif

   open(unit=stderrUnit,file=trim(errorfile))
   write(stderrUnit,*) 'MPAS critical error'
   write(stderrUnit,*) trim(mesg)
   close(stderrUnit)

   ! stop the code
#ifdef _MPI
   call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#else
   stop
#endif

end subroutine mpas_dmpar_global_abort!}}}

