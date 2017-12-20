
#include "mpi.h"
#ifndef ZOLTANINTERFACEHPP
#define ZOLTANINTERFACEHPP


#ifdef __cplusplus
extern "C" {
#endif
void sort_graph(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt);
#ifdef __cplusplus
}
#endif


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
     double *zcoord, int *coord_dimension,
     int *result_parts,
     int *partmethod,
     int *mapmethod);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
 void zoltan_map_problem(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     double *xcoord,
     double *ycoord,
     double *zcoord, int *coord_dimension,
     int *result_parts,
     int *partmethod,
     int *mapmethod);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
 void zoltan_mapping_problem(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     double *xcoord,
     double *ycoord,
     double *zcoord, int *coord_dimension,
     int *result_parts,
     int *partmethod,
     int *mapmethod);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
 void zoltan2_print_metrics(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     int *result_parts);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
 void zoltan2_print_metrics2(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     int *result_parts);
#ifdef __cplusplus
}
#endif


#endif
