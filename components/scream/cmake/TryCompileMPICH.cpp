#include <mpi.h>

int main (int, char**)
{
#ifndef MPICH_NAME
#error "Mpi not provided by OpenMPI"
#endif
  return 0;
}

