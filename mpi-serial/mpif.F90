#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

Module mpi
implicit none
! MPI_ADDRESS_KIND: typekind 4 should be 4 byte integer,
!                   and typekind 8 be 8-byte integer, etc.
#ifdef SIZEOF_LONG
        INTEGER, PARAMETER, PUBLIC :: MPI_ADDRESS_KIND=SIZEOF_LONG
#else
        INTEGER, PARAMETER, PUBLIC :: MPI_ADDRESS_KIND=8
#endif


        include "mpif.h"
end Module mpi
