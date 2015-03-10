program main 
  use gptl

  implicit none

#include <mpif.h>

#ifdef THREADED_OMP
  integer, external :: omp_get_max_threads
#endif

  double precision, external :: sub
  double precision result

  integer :: iam
  integer :: nthreads = 1 ! number of threads (default 1)
  integer :: nproc = 1
  integer iter
  integer code
  integer c
  integer :: comm = 0
  integer ierr
  integer ret

#ifdef HAVE_PAPI
! Turn abort_on_error off just long enough to check PAPI-based options
  ret = gptlsetoption (gptlabort_on_error, 0)
  if (gptlevent_name_to_code ('PAPI_FP_OPS', code) == 0) then
    ret = gptlsetoption (code, 1)
  end if
  ret = gptlsetoption (gptl_ci, 1)
  ret = gptlsetoption (gptlabort_on_error, 1)
#endif

  ret = gptlsetoption (gptlabort_on_error, 1)
  ret = gptlsetoption (gptloverhead, 1)
  ret = gptlsetoption (gptlnarrowprint, 1)

  call mpi_init (ierr)
  comm = MPI_COMM_WORLD

#ifndef ENABLE_PMPI
  ret = gptlinitialize ()
  ret = gptlstart ("total")
#endif
	 
  call mpi_comm_rank (MPI_COMM_WORLD, iam, ierr)
  call mpi_comm_size (MPI_COMM_WORLD, nproc, ierr)

  if (iam == 0) then
    write (6,*) "Purpose: test behavior of summary stats"
    write (6,*) "Include OpenMP if enabled"
  end if

#ifdef THREADED_OMP
  nthreads = omp_get_max_threads ()
#endif

!$OMP PARALLEL DO PRIVATE (RESULT)
  do iter=1,nthreads
    result = sub (iter, iam)
  end do

#ifndef ENABLE_PMPI
  ret = gptlstop ("total")
  ret = gptlpr (iam)
#endif
  ret = gptlpr_summary (comm)
  if (ret /= 0) then
    write(6,*)'summary.F90: error from gptlpr_summary'
    stop 1
  end if
  ret = gptlpr_summary_file (comm, "timing.summary.duplicate")
  if (ret /= 0) then
    write(6,*)'summary.F90: error from gptlpr_summary_file'
    stop 1
  end if

  call mpi_finalize (ret)

  if (gptlfinalize () < 0) stop 1
  stop 0
end program main


double precision function sub (iter, iam)
  use gptl
  implicit none
  
  integer, intent (in) :: iter
  integer, intent (in) :: iam

  integer (8) :: looplen
  integer (8) :: i
  integer :: ret
  double precision sum

  looplen = iam*iter*10000
  ret = gptlstart ("sub")

  ret = gptlstart ("sleep")
  ret = gptlstop ("sleep")

  ret = gptlstart ("work")
  sum = 0.
  ret = gptlstart ("add")
  do i=0,looplen-1
    sum = sum + i
  end do
  ret = gptlstop ("add")

  ret = gptlstart ("madd")
  do i=0,looplen-1
    sum = sum + i*1.1
  end do
  ret = gptlstop ("madd")
  
  ret = gptlstart ("div")
  do i=0,looplen-1
    sum = sum / 1.1
  end do
  ret = gptlstop ("div")
  ret = gptlstop ("work")
  ret = gptlstop ("sub")
  
  sub = sum
  return 
end function sub
