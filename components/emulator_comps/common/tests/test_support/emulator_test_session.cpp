//==============================================================================
// Emulator Test Session
//
// Custom Catch2 main() with MPI initialization/finalization.
// All tests run across MPI ranks, with results aggregated.
// Note: CATCH_CONFIG_RUNNER is defined by CMake via target_compile_definitions
//==============================================================================

#include <catch2/catch.hpp>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  // Initialize MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Create Catch2 session
  Catch::Session session;

  // Parse command line arguments
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) {
    MPI_Finalize();
    return returnCode;
  }

  // Print header on rank 0
  if (rank == 0) {
    std::cout << "========================================" << std::endl;
    std::cout << "Emulator Component Tests" << std::endl;
    std::cout << "MPI ranks: " << nprocs << std::endl;
    std::cout << "========================================" << std::endl;
  }

  // Synchronize before running tests
  MPI_Barrier(MPI_COMM_WORLD);

  // Run tests
  int result = session.run();

  // Synchronize after tests
  MPI_Barrier(MPI_COMM_WORLD);

  // Aggregate results across all ranks (if any rank failed, report failure)
  int global_result;
  MPI_Allreduce(&result, &global_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Print summary on rank 0
  if (rank == 0) {
    std::cout << "========================================" << std::endl;
    if (global_result == 0) {
      std::cout << "All tests passed on all ranks!" << std::endl;
    } else {
      std::cout << "Tests failed on one or more ranks." << std::endl;
    }
    std::cout << "========================================" << std::endl;
  }

  MPI_Finalize();

  return global_result;
}
