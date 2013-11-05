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

   use mpas_dmpar_types
   use mpas_grid_types
   use mpas_sort
   use mpas_hash
   use mpas_io_units

#ifdef _MPI
include 'mpif.h'
   integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER

#ifdef SINGLE_PRECISION
   integer, parameter :: MPI_REALKIND = MPI_REAL
#else
   integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
#endif
#endif

   integer, parameter :: IO_NODE = 0
   integer, parameter :: BUFSIZE = 6000

   interface mpas_dmpar_alltoall_field
      module procedure mpas_dmpar_alltoall_field1d_integer
      module procedure mpas_dmpar_alltoall_field2d_integer
      module procedure mpas_dmpar_alltoall_field1d_real
      module procedure mpas_dmpar_alltoall_field2d_real
      module procedure mpas_dmpar_alltoall_field3d_real
      module procedure mpas_dmpar_alltoall_field4d_real
      module procedure mpas_dmpar_alltoall_field5d_real
   end interface

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

   private :: mpas_dmpar_copy_field1d_integer
   private :: mpas_dmpar_copy_field2d_integer
   private :: mpas_dmpar_copy_field3d_integer
   private :: mpas_dmpar_copy_field1d_real
   private :: mpas_dmpar_copy_field2d_real
   private :: mpas_dmpar_copy_field3d_real
   private :: mpas_dmpar_copy_field4d_real
   private :: mpas_dmpar_copy_field5d_real

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

      if (present(mpi_comm)) then
         dminfo % comm = mpi_comm
         dminfo % using_external_comm = .true.
      else
         call MPI_Init(mpi_ierr)
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

#ifndef MPAS_CESM
      call open_streams(dminfo % my_proc_id)
#endif

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
!  routine mpas_dmpar_global_abort
!
!> \brief MPAS dmpar global abort routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine aborts MPI. A call to it kills the model through the use of MPI_Abort on the world communicator, and outputs a message.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_global_abort(mesg)!{{{

      implicit none

      character (len=*), intent(in) :: mesg !< Input: Abort message

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      write(stderrUnit,*) trim(mesg)
      call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#endif

      write(stderrUnit,*) trim(mesg)
      stop

   end subroutine mpas_dmpar_global_abort!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_int
!
!> \brief MPAS dmpar broadcast integer routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts an integer to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_int(dminfo, i)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(inout) :: i !< Input/Output: Integer to broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(i, 1, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_int!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_ints
!
!> \brief MPAS dmpar broadcast integers routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts an array of integers to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_ints(dminfo, n, iarray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      integer, dimension(n), intent(inout) :: iarray !< Input/Output: Array of integers

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(iarray, n, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_ints!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_real
!
!> \brief MPAS dmpar broadcast real routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts a real to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_real(dminfo, r)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      real (kind=RKIND), intent(inout) :: r !< Input/Output: Real to be broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(r, 1, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_real!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_reals
!
!> \brief MPAS dmpar broadcast reals routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts an array of reals to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_reals(dminfo, n, rarray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      real (kind=RKIND), dimension(n), intent(inout) :: rarray !< Input/Output: Array of reals to be broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(rarray, n, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_reals!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_double
!
!> \brief MPAS dmpar broadcast double routine.
!> \author Michael Duda
!> \date   11/04/13
!> \details
!>  This routine broadcasts a double to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_double(dminfo, r)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      double precision, intent(inout) :: r !< Input/Output: Double to be broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(r, 1, MPI_DOUBLE_PRECISION, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_double!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_doubles
!
!> \brief MPAS dmpar broadcast doubles routine.
!> \author Michael Duda
!> \date   11/04/13
!> \details
!>  This routine broadcasts an array of doubles to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_doubles(dminfo, n, rarray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      integer, intent(in) :: n !< Input: Length of array
      double precision, dimension(n), intent(inout) :: rarray !< Input/Output: Array of doubles to be broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(rarray, n, MPI_DOUBLE_PRECISION, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_doubles!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_logical
!
!> \brief MPAS dmpar broadcast logical routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts a logical to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_logical(dminfo, l)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      logical, intent(inout) :: l !< Input/Output: Logical to be broadcast

#ifdef _MPI
      integer :: mpi_ierr
      integer :: itemp

      if (dminfo % my_proc_id == IO_NODE) then
         if (l) then
            itemp = 1
         else
            itemp = 0
         end if
      end if

      call MPI_Bcast(itemp, 1, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)

      if (itemp == 1) then
         l = .true.
      else
         l = .false.
      end if
#endif

   end subroutine mpas_dmpar_bcast_logical!}}}

!-----------------------------------------------------------------------
!  routine mpas_dmpar_bcast_char
!
!> \brief MPAS dmpar broadcast character routine.
!> \author Michael Duda
!> \date   03/26/13
!> \details
!>  This routine broadcasts a character to all processors in the communicator.
!
!-----------------------------------------------------------------------
   subroutine mpas_dmpar_bcast_char(dminfo, c)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      character (len=*), intent(inout) :: c !< Input/Output: Character to be broadcast

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(c, len(c), MPI_CHARACTER, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_char!}}}

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

#ifdef _MPI
      call MPI_Allreduce(i, isum, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      isum = i
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rsum, 1, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      rsum = r
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(i, imin, 1, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      imin = i
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rmin, 1, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      rmin = r
#endif

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
      
      integer :: mpi_ierr 
      
#ifdef _MPI
      call MPI_Allreduce(i, imax, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      imax = i
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rmax, 1, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      rmax = r
#endif

   end subroutine mpas_dmpar_max_real!}}}

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
      
      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

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

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_max_real_array!}}}

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
      integer :: mpi_ierr
      
      call MPI_Scatterv(inlist, counts, displs, MPI_INTEGERKIND, outlist, noutlist, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
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

      local_start = nint(real(dminfo % my_proc_id) * real(global_end - global_start + 1) / real(dminfo % nprocs)) + 1
      local_end   = nint(real(dminfo % my_proc_id + 1) * real(global_end - global_start + 1) / real(dminfo % nprocs)) 

   end subroutine mpas_dmpar_get_index_range!}}}

   subroutine mpas_dmpar_compute_index_range(dminfo, &!{{{
                                        local_start, local_end, &
                                        global_start, global_end)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: local_start, local_end
      integer, intent(inout) :: global_start, global_end

      integer :: n
      integer :: mpi_ierr

      n = local_end - local_start + 1

      if (dminfo % my_proc_id == 0) then
         global_start = 1
         global_end = global_start + n - 1
         
#ifdef _MPI
      else if (dminfo % my_proc_id == dminfo % nprocs - 1) then
         call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         global_end = global_start + n - 1

      else
         call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         global_end = global_start + n
         call MPI_Send(global_end, 1, MPI_INTEGERKIND, dminfo % my_proc_id + 1, 0, dminfo % comm, mpi_ierr)
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
      integer :: totalSize, nMesgRecv, nMesgSend, recvNeighbor, sendNeighbor, currentProc, offset
      integer :: totalSent, totalRecv
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
      integer :: nUniqueNeededList
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

      ! 
      ! Determine total number of owned blocks on this task
      ! 
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
      deallocate(ownedList)
      deallocate(ownedListIndex)
      deallocate(ownedBlock)
      deallocate(neededList)
      deallocate(neededListIndex)
      deallocate(neededBlock)
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
   subroutine mpas_dmpar_build_comm_lists(sendExchList, recvExchList, haloLayers, dimsizes, sendCommList, recvCommList)

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
      type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
      type (mpas_communication_list), pointer :: commListPtr, commListPtr2
      logical :: comm_list_found
      integer :: nAdded, bufferOffset
      integer :: iHalo
      integer :: nHaloLayers
      integer :: iDimen
      integer :: dimSizeProduct  ! the product of the size of all dimensions


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

   !--------------------------------------------------------------------
   end subroutine mpas_dmpar_build_comm_lists



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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k, l
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k, l, m
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k, l
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 4
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k, l, m
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 5
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

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
                          fieldCursor % array(m, l, k, j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)&
                                                                                 *fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) *fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4)&
                                                                               + (j-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) &
                                                                               + (k-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                               + (l-1)*fieldCursor % dimSizes(1) + m + bufferOffset)
                        end do
                      end do
                    end do
                  end do
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) * fieldCursor % dimSizes(3) * fieldCursor % dimSizes(4))
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

     integer :: i

     allocate(exchList)
     allocate(exchList % halos(nHalos))
     do i = 1, nHalos
       nullify(exchList % halos(i) % exchList)
     end do
     nullify(exchList % next)
     nullify(exchList % prev)
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
     integer :: i

     nHalos = size(exchList % halos)

     do i = 1, nHalos
       call mpas_dmpar_destroy_exchange_list(exchList % halos(i) % exchList)
     end do

     deallocate(exchList % halos)
     deallocate(exchList)
     nullify(exchList)
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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


       if(associated(field % next)) then
         fieldCursor => field
         do while(associated(fieldCursor))
           fieldCursor % array(:) = field % array(:)
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
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

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field5d_real!}}}

end module mpas_dmpar
