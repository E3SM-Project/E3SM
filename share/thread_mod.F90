#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module thread_mod
  implicit none
  private

  integer, public :: NThreads, nthreads_accel

#ifdef _OPENMP
  interface omp_get_thread_num
     integer function omp_get_thread_num()
     end function omp_get_thread_num
  end interface

  interface omp_in_parallel
     function omp_in_parallel() result(ans)
       logical ans
     end function omp_in_parallel
  end interface

  interface omp_set_num_threads
     subroutine omp_set_num_threads(NThreads)
       integer NThreads
     end subroutine omp_set_num_threads
  end interface

  interface omp_get_num_threads
     integer function omp_get_num_threads()
       integer NThreads
     end function omp_get_num_threads
  end interface

  interface omp_get_max_threads
     integer function omp_get_max_threads()
     end function omp_get_max_threads
  end interface
#endif
  public :: omp_get_thread_num
  public :: omp_in_parallel
  public :: omp_set_num_threads
  public :: omp_get_max_threads
  public :: omp_get_num_threads
#ifndef _OPENMP
contains

  function omp_get_thread_num() result(ithr)
    integer ithr
    ithr=0
  end function omp_get_thread_num

  function omp_get_num_threads() result(ithr)
    integer ithr
    ithr=1
  end function omp_get_num_threads

  function omp_in_parallel() result(ans)
    logical ans
    ans=.FALSE.
  end function omp_in_parallel

  subroutine omp_set_num_threads(NThreads)
    integer Nthreads
    NThreads=1
  end subroutine omp_set_num_threads

  integer function omp_get_max_threads()
    omp_get_max_threads=1
  end function omp_get_max_threads
#endif

end module thread_mod
