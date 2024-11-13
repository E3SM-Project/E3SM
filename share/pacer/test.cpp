#include <iostream>
#include <mpi.h>
#include "Pacer.h"

int main(int argc, char **argv){

    int err;
    int myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    Pacer::initialize(MPI_COMM_WORLD);

    Pacer::setPrefix("Omega");

    Pacer::start("run_loop");

    float tmp = 1;

    for (int i = 1; i <= 1000; i++){
        tmp *= i;
    }

    Pacer::stop("run_loop");


    if (myrank == 0)
        Pacer::print("omega");

    Pacer::finalize();
}

