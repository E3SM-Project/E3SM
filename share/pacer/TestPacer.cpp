////===-- TestPacer.cpp - Simple test for Pacer --*- C++ -*-===//
//
// \file
// \brief Simple example illustrating Pacer API usage.
//
// This test exercises basic timer functionality
// with the Pacer API.
//
// This test program should create two files:
// test.timing.0 and test.summary
// It is also expected to issue couple of warnings
// to illustrate likely scenarios where a timer is
// not started/stopped properly.
//
////===-----------------------------------------------------===//

#include <iostream>
#include <mpi.h>
#include "Pacer.h"

int main(int argc, char **argv){

    int err;
    int myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    Pacer::initialize(MPI_COMM_WORLD);
    // Second argument is optional (default is Pacer::PACER_STANDALONE)
    // Pacer::initialize(MPI_COMM_WORLD, Pacer::PACER_STANDALONE);

    Pacer::setPrefix("Omega:");

    Pacer::start("run_loop");

    float tmp = 1;

    for (int i = 1; i <= 10000; i++){
        tmp *= i;
    }

    Pacer::stop("run_loop");

    Pacer::unsetPrefix();

    // illustrating situation where attempt to stop timer before starting
    // will print a warning
    Pacer::stop("final");

    // as well as dangling timer which is never stopped explicitly
    // will print a warning at the end during finalize
    Pacer::start("final");


    // print: First argument is the prefix used for timing output file names
    // Second argument (optional) controls if timing output should be from all ranks
    // default is false, only rank 0 writes timing output
    Pacer::print("test");
    // Pacer::print("test", true);

    Pacer::finalize();
}

