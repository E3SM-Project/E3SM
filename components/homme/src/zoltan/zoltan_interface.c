
#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#if HAVE_TRILINOS
#if TRILINOS_HAVE_ZOLTAN2
#include "zoltan_cppinterface.hpp"
#endif
#endif
#include "stdio.h"

//#define VISUALIZEOUTPUT

void zoltanpart_(int *nelem, int *xadj,int *adjncy,double *adjwgt,double *vwgt, int *nparts, MPI_Fint *comm,
    double *xcoord, double *ycoord, double *zcoord, int *result_parts, int *partmethod) {
  MPI_Comm c_comm = MPI_Comm_f2c(*comm);


#ifdef VISUALIZEINPUT
  int mype2, size2;
  MPI_Comm_rank(c_comm, &mype2);
  MPI_Comm_size(c_comm, &size2);
  if (mype2 == 0) {
    int i = 0, j =0;
    FILE *f2 = fopen("preplot.gnuplot", "w");

    FILE *coord_files = fopen("preplotcoords.txt", "w");

    fprintf(f2,"plot \"preplotcoords.txt\"\n");
    for (j = 0; j < *nelem; ++j){
      fprintf(coord_files,"%lf %lf %lf\n", xcoord[j], ycoord[j], zcoord[j]);
    }

    fclose(coord_files);



    for (i = 0; i < *nelem; ++i){
      for (j = xadj[i]; j < xadj[i+1] ; ++j){
        int n = adjncy[j];
        fprintf(f2,"set arrow from %lf,%lf,%lf to %lf,%lf,%lf nohead lt -1 lw 1.2\n", xcoord[i], ycoord[i], zcoord[i],
                    xcoord[n], ycoord[n], zcoord[n]);
      }
    }

    fprintf(f2,"pause-1\n");
    fclose(f2);


  }
#endif

#if HAVE_TRILINOS
#if TRILINOS_HAVE_ZOLTAN2
  zoltan_partition_problem(nelem, xadj,adjncy,adjwgt,vwgt, nparts, c_comm,
      xcoord, ycoord, zcoord, result_parts, partmethod);
#else
  int mype2, size2;
    MPI_Comm_rank(c_comm, &mype2);
    MPI_Comm_size(c_comm, &size2);
    if (mype2 == 0) {

      printf("Zoltan cannot be used, Trilinos is not compiled with Zoltan2.");
    }
    exit(1);
#endif
#else
  int mype2, size2;
  MPI_Comm_rank(c_comm, &mype2);
  MPI_Comm_size(c_comm, &size2);
  if (mype2 == 0) {

    printf("Zoltan cannot be used, HOMME is not compiled with Trilinos.");
  }
  exit(1);
#endif

#ifdef VISUALIZEOUTPUT
  int mype, size;
  MPI_Comm_rank(c_comm, &mype);
  MPI_Comm_size(c_comm, &size);
  if (mype == 0) {
    int i = 0, j =0;
    FILE *f2 = fopen("plot.gnuplot", "w");

    FILE **coord_files = (FILE **) malloc(sizeof(FILE*) * size);
    for (i = 0; i< size; ++i){
      char str[20];
      sprintf(str, "coords%d.txt", i);
      coord_files[i] = fopen(str, "w");
      if (i == 0){
        fprintf(f2,"splot \"%s\"\n",  str);
      }
      else {
        fprintf(f2,"replot \"%s\"\n",  str);
      }
    }

    for (j = 0; j < *nelem; ++j){
      int findex = result_parts[j];
      fprintf(coord_files[findex - 1],"%lf %lf %lf\n", xcoord[j], ycoord[j], zcoord[j]);
    }
    for (int i = 0; i< size; ++i){
      fclose(coord_files[i]);
    }


    fprintf(f2,"pause-1\n");
    for (i = 0; i < *nelem; ++i){
      for (j = xadj[i]; j < xadj[i+1] ; ++j){
        int n = adjncy[j];
        fprintf(f2,"set arrow from %lf,%lf,%lf to %lf,%lf,%lf nohead lt -1 lw 1.2\n", xcoord[i], ycoord[i], zcoord[i],
                    xcoord[n], ycoord[n], zcoord[n]);
      }
    }

    fprintf(f2,"pause-1\n");
    fclose(f2);
    free(coord_files);

  }
#endif

}
