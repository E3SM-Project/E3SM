#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include "share/scream_session.hpp"
#include <mpi.h>

int main (int argc, char **argv) {

  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Initialize scream;
  scream::initialize_scream_session(argc, argv);

  // Run tests
  int result = Catch::Session().run(argc, argv);
  
  // Finalize scream
  scream::finalize_scream_session();

  // Finalize MPI
  MPI_Finalize();

  // Return test result
  return result;
}
