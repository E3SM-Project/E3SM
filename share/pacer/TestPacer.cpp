// This test exercises basic timer functionality
// with the Pacer API.
//
// This test program should create two files:
// test_pacer.timing and test_pacer.summary
// It is also expected to issue a warning about 
// the "final" timer still being open.

#include <iostream>
#include <mpi.h>
#include "Pacer.h"

int main(int argc, char **argv){

    int err;
    int myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    Pacer::initialize(MPI_COMM_WORLD);

    Pacer::setPrefix("Omega:");

    Pacer::start("run_loop");

    float tmp = 1;

    for (int i = 1; i <= 10000; i++){
        tmp *= i;
    }

    Pacer::stop("run_loop");

    Pacer::unsetPrefix();
    Pacer::start("final");

    if (myrank == 0)
        Pacer::print("test_pacer");

    Pacer::finalize();
}

