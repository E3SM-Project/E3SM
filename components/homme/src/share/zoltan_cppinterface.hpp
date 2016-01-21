
#include "mpi.h"
#ifndef ZOLTANINTERFACEHPP
#define ZOLTANINTERFACEHPP
#ifdef __cplusplus
extern "C" {
#endif
 void zoltan_partition_problem(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     double *xcoord,
     double *ycoord,
     double *zcoord,
     int *result_parts,
     int *partmethod);
#ifdef __cplusplus
}
#endif

#endif
