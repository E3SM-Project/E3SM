#include <mpi.h>

int main (int, char**)
{
#ifdef MPICH_VERSION
#error "MPICH found"
#endif
  return 0;
}

