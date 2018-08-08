
#define CATCH_CONFIG_RUNNER

#include "catch/catch.hpp"

#include <Kokkos_Core.hpp>

#include <Hommexx_Session.hpp>
#include <Context.hpp>
#include <mpi/Comm.hpp>

#include <mpi.h>

int main(int argc, char **argv) {

  // Initialize mpi
  MPI_Init(&argc,&argv);

  Homme::initialize_hommexx_session();
  Homme::Context::singleton().get_comm().reset_mpi_comm(MPI_COMM_WORLD);

  int result = Catch::Session().run(argc, argv);

  Homme::finalize_hommexx_session();

  // Finalize mpi
  MPI_Finalize();

  return result;
}
