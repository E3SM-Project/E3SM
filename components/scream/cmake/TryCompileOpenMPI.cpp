#include <mpi.h>

int main (int, char**)
{
#ifdef OMPI_MAJOR_VERSION
#error "OpenMPI found"
#endif
  return 0;
}
