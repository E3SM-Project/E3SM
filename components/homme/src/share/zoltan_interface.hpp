
#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif


#ifdef HAVE_ZOLTAN
#include "Zoltan2_TaskMapping.hpp"
//#include "Teuchos_Comm.hpp"
//#include "Teuchos_CommHelpers.hpp"
//#include "Teuchos_RCPDecl.hpp"
#else
#endif
#include "mpi.h"
#include <stdio.h>

//#include <iostream> 
#include "stdio.h"
#define VISUALIZE
extern "C" void ZOLTANPART_(int *nelem, int *xadj,int *adjncy,int *adjwgt,int *vwgt, int *nparts, MPI_Fint *comm,
    double *xcoord, double *ycoord, double *zcoord) {

}
extern "C"
{
void zoltanpart_(int *nelem, int *xadj,int *adjncy,int *adjwgt,int *vwgt, int *nparts, MPI_Fint *comm,
    double *xcoord, double *ycoord, double *zcoord) {

  int mype;
  MPI_Comm c_comm = MPI_Comm_f2c(*comm);
  MPI_Comm_rank(c_comm, &mype);

#ifdef VISUALIZE
  if (mype == 0) {
#ifdef HAVE_ZOLTAN
    printf("\n\n\n\nZOLTAN IS DEFINED!!!!\n\n\n\n\n");
#else
    printf("\n\n\n\nZOLTAN ISnot DEFINED!!!!\n\n\n\n\n");
#endif
    FILE *f = fopen("coords.txt", "w");
    FILE *f2 = fopen("plot.gnuplot", "w");
    fprintf(f2,"splot \"coords.txt\"\n");

    printf("ZOLTAN nelem:%d\n", *nelem);
    printf("ZOLTAN nedge:%d\n", xadj[*nelem]);
    printf("ZOLTAN nparts:%d\n", *nparts);
    int i = 0, j =0;
    for (i = 0; i < *nelem; ++i){
      if (i == 1)
      printf("i:%d adjs:%d adje:%d %lf %lf %lf\n\t",i, xadj[i], xadj[i+1], xcoord[i], ycoord[i], zcoord[i]);

      for (int j = xadj[i]; j < xadj[i+1] ; ++j){
        //printf("n:%d ", adjncy[j]);
        int n = adjncy[j];
        if (i % 16 == 0 || 1){
          //printf("n:%d ", adjncy[j]);
          fprintf(f2,"set arrow from %lf,%lf,%lf to %lf,%lf,%lf nohead lt -1 lw 1.2\n", xcoord[i], ycoord[i], zcoord[i],
            xcoord[n], ycoord[n], zcoord[n]);
        }
      }
      //printf("\n");

      fprintf(f2,"set label \"%d\" at %lf,%lf,%lf\n", i, xcoord[i], ycoord[i], zcoord[i]);
      fprintf(f,"%lf %lf %lf\n", xcoord[i], ycoord[i], zcoord[i]);
    }


    fprintf(f2,"pause-1\n");
    fclose(f2);
    fclose(f);
  }
#endif
  //std::cout << "IN C++ FILE" << std::endl;
/*
  Zoltan2::coordinateTaskMapperInterface<int, float, float>(
	tcomm,
        3,
        numNodes,
        coords,

        3,
        numNodes,
        taskcoords,

        task_communication_xadj_,
        task_communication_adj_,

        );
*/

}
}


//#endif
