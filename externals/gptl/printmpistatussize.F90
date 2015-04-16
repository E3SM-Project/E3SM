program printmpistatussize
#ifdef HAVE_MPI
  use mpi
#endif
  implicit none
#ifdef HAVE_MPI
  write(6,*) 'MPI_STATUS_SIZE=', MPI_STATUS_SIZE
#else
  write(6,*) 'Need set HAVE_MPI=yes to get MPI_STATUS_SIZE'
#endif
end program printmpistatussize
