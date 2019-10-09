program mpicheck
  include "mpif.h"
  integer :: fh, ierr

  call mpi_file_open(mpi_comm_world, 'stupid.file',MPI_MODE_RDWR,MPI_INFO_NULL,fh,ierr)
end program
