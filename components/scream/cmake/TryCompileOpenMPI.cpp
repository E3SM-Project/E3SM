#include <mpi.h>

int main (int, char**)
{
#ifndef OMPI_MAJOR_VERSION
#error "Mpi not provided by OpenMPI"
#endif
  return 0;
}
