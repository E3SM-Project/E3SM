module thread_mod

#ifdef _OPENMP
  use omp_lib, only: omp_get_thread_num, &
       omp_in_parallel, &
       omp_set_num_threads, &
       omp_get_max_threads, &
       omp_get_num_threads, &
       omp_get_nested
#endif
  use cam_logfile,            only: iulog
  use spmd_utils,             only: masterproc

  implicit none
  private

  integer, public :: max_num_threads=1 ! maximum number of OpenMP threads
  integer, public :: horz_num_threads, vert_num_threads, tracer_num_threads

  public :: omp_get_thread_num
  public :: omp_in_parallel
  public :: omp_set_num_threads
  public :: omp_get_max_threads
  public :: omp_get_num_threads
  public :: omp_get_nested
  public :: initomp
contains

#ifndef _OPENMP
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

  integer function omp_get_nested()
    omp_get_nested=0
  end function omp_get_nested

  subroutine initomp
    max_num_threads = 1
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: openmp not activated"
    end if
  end subroutine initomp

#else
  subroutine initomp
    max_num_threads = omp_get_num_threads()
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: number of OpenMP threads = ", max_num_threads
    end if
  end subroutine initomp
#endif

end module thread_mod
