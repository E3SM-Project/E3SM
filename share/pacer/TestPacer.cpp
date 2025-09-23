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
    Pacer::setTimingLevel(1);

    Pacer::start("run_loop", 1);

    float tmp = 1;

    for (int i = 1; i <= 10000; i++){
        tmp *= i;
    }

    Pacer::stop("run_loop", 1);
    
    Pacer::start("should_not_appear_1", 2);
    Pacer::stop("should_not_appear_1", 2);
    
    Pacer::disableTiming();
      Pacer::start("should_not_appear_2", 0);

      Pacer::start("should_not_appear_3", 0);
      Pacer::stop("should_not_appear_3", 0);

      Pacer::stop("should_not_appear_2", 0);
    Pacer::enableTiming();
    
    Pacer::start("parent", 1);
    Pacer::addParentPrefix();
    
    Pacer::start("child1", 1);
    Pacer::stop("child1", 1);
    
    Pacer::start("child2", 0);
      Pacer::addParentPrefix();
      Pacer::start("grandchild1", 1);
      Pacer::stop("grandchild1", 1);
      Pacer::removeParentPrefix();
    Pacer::stop("child2", 0);

    Pacer::removeParentPrefix();
    Pacer::stop("parent", 1);

    Pacer::timingBarrier("barrier_should_not_appear_1", 0, MPI_COMM_WORLD);
    Pacer::enableTimingBarriers();
    Pacer::timingBarrier("barrier_should_not_appear_2", 2, MPI_COMM_WORLD);
    Pacer::timingBarrier("timing_barrier", 0, MPI_COMM_WORLD);


    Pacer::unsetPrefix();

    // illustrating situation where attempt to stop timer before starting
    // will print a warning
    Pacer::stop("final", 0);

    // as well as dangling timer which is never stopped explicitly
    // will print a warning at the end during finalize
    Pacer::start("final", 0);


    // print: First argument is the prefix used for timing output file names
    // Second argument (optional) controls if timing output should be from all ranks
    // default is false, only rank 0 writes timing output
    Pacer::print("pacer_test");
    // Pacer::print("test", true);

    Pacer::finalize();
}

