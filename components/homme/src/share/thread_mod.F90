#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module thread_mod

#ifdef _OPENMP
  use omp_lib, only: omp_get_thread_num, &
       omp_in_parallel, &
       omp_set_num_threads, &
       omp_get_max_threads, &
       omp_get_num_threads, &
       omp_get_nested
#endif
#ifdef CAM
  use cam_logfile,            only: iulog
  use spmd_utils,             only: masterproc
#endif

  implicit none
  private

#ifdef CAM
  integer, public,  TARGET :: max_num_threads ! maximum number of OpenMP threads
  integer, public          :: tracer_num_threads
  integer, public,  TARGET :: horz_num_threads , vert_num_threads

  integer, public, pointer :: NThreads   ! total number of threads
                                         ! standalone HOMME: from namelist
                                         ! in CAM: set by driver
  integer, public, pointer :: hthreads   ! computed based on nthreads, vthreads,nelemd
  integer, public, pointer :: vthreads   ! not used unless set in namelist
#else
  integer, public :: NThreads   ! total number of threads
                                ! standalone HOMME: from namelist
                                ! in CAM: set by driver
  integer, public :: hthreads   ! computed based on nthreads, vthreads,nelemd
  integer, public :: vthreads = 1   ! not used unless set in namelist
#endif
  public :: omp_get_thread_num
  public :: omp_in_parallel
  public :: omp_set_num_threads
  public :: omp_get_max_threads
  public :: omp_get_num_threads
  public :: omp_get_nested
#ifdef CAM
  public :: initomp
#endif

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

  integer function omp_get_nested()
    omp_get_nested=0
  end function omp_get_nested

#ifdef CAM
  subroutine initomp
    max_num_threads = 1
    NThreads=>max_num_threads
    hthreads=>horz_num_threads
    vthreads => vert_num_threads
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: openmp not activated"
    end if
  end subroutine initomp
#endif

#else
#ifdef CAM
contains

  subroutine initomp
    !$OMP PARALLEL
    max_num_threads = omp_get_num_threads()
    !$OMP END PARALLEL
    NThreads=>max_num_threads
    hthreads=>horz_num_threads
    vthreads => vert_num_threads
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: number of OpenMP threads = ", max_num_threads
    end if
  end subroutine initomp
#endif
#endif

end module thread_mod
