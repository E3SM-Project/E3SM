! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_threading
!
!> \brief MPAS Threading Support
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This module will provide interfaces to support functions / routines for OpenMP threading.
!
!-----------------------------------------------------------------------

module mpas_threading

   use mpas_kind_types

   implicit none
   private

   public :: mpas_threading_get_num_threads, mpas_threading_set_num_threads, mpas_threading_in_parallel
   public :: mpas_threading_get_thread_num, mpas_threading_barrier, mpas_threading_get_max_threads 
   public :: mpas_threading_get_thread_limit

   contains

!-----------------------------------------------------------------------
!  function mpas_threading_get_num_threads
!
!> \brief MPAS Threading number of threads function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns the number of threads currently available.
!
!-----------------------------------------------------------------------
   function mpas_threading_get_num_threads() result(numThreads)!{{{
      integer :: numThreads
      integer :: omp_get_num_threads

      numThreads = 1

#ifdef MPAS_OPENMP
      numThreads = omp_get_num_threads()
#endif

   end function mpas_threading_get_num_threads!}}}


!-----------------------------------------------------------------------
!  routine mpas_threading_set_num_threads
!
!> \brief MPAS Threading set number of threads routine
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This routine sets the number of threads for the next parallel region.
!
!-----------------------------------------------------------------------
   subroutine mpas_threading_set_num_threads(numThreads)!{{{
      integer, intent(in) :: numThreads

#ifdef MPAS_OPENMP
      call omp_set_num_threads(numThreads)
#endif

   end subroutine mpas_threading_set_num_threads!}}}

!-----------------------------------------------------------------------
!  function mpas_threading_in_parallel
!
!> \brief MPAS Threading in parallel function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns a logical where true means it was called within a
!>  parallel region, and false means it was not.
!
!-----------------------------------------------------------------------
   function mpas_threading_in_parallel() result(parallelRegion)!{{{
      logical :: parallelRegion
      logical :: omp_in_parallel

      parallelRegion = .false.

#ifdef MPAS_OPENMP
      parallelRegion = omp_in_parallel()
#endif

   end function mpas_threading_in_parallel!}}}

!-----------------------------------------------------------------------
!  function mpas_threading_get_thread_num
!
!> \brief MPAS Threading get thread number function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns current thread's number
!
!-----------------------------------------------------------------------
   function mpas_threading_get_thread_num() result(threadNum)!{{{
      integer :: threadNum
      integer :: omp_get_thread_num

      threadNum = 0

#ifdef MPAS_OPENMP
      threadNum = omp_get_thread_num()
#endif

   end function mpas_threading_get_thread_num!}}}

!-----------------------------------------------------------------------
!  routine mpas_threading_barrier
!
!> \brief MPAS Threading barrier routine
!> \author Doug Jacobsen
!> \date   10/15/2015
!> \details
!>  This routine implements an OpenMP barrier to synchronize all threads.
!
!-----------------------------------------------------------------------
   subroutine mpas_threading_barrier()!{{{

#ifdef MPAS_OPENMP
      !$omp barrier
#endif

   end subroutine mpas_threading_barrier!}}}

!-----------------------------------------------------------------------
!  function mpas_threading_get_max_threads
!
!> \brief MPAS Threading maximum number of threads function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns maximum number of threads a single MPI process can use.
!
!-----------------------------------------------------------------------
   function mpas_threading_get_max_threads() result(maxThreads)!{{{
      integer :: maxThreads
      integer :: omp_get_max_threads

      maxThreads = 1

#ifdef MPAS_OPENMP
      maxThreads = omp_get_max_threads()
#endif

   end function mpas_threading_get_max_threads!}}}

!-----------------------------------------------------------------------
!  function mpas_threading_get_thread_limit
!
!> \brief MPAS Threading thread limit function
!> \author Doug Jacobsen
!> \date   09/09/2015
!> \details
!>  This function returns limit on the total number of threads.
!
!-----------------------------------------------------------------------
   function mpas_threading_get_thread_limit() result(threadLimit)!{{{
      integer :: threadLimit
      integer :: omp_get_thread_limit

      threadLimit = 1

#ifdef MPAS_OPENMP
      threadLimit = omp_get_thread_limit()
#endif

   end function mpas_threading_get_thread_limit!}}}

end module mpas_threading
