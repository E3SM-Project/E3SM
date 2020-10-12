
#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include <Kokkos_Core.hpp>

#include <Hommexx_Session.hpp>

#include <Context.hpp>
#include <mpi/Comm.hpp>

#include <mpi.h>

int main(int argc, char **argv) {

  // Initialize mpi
  MPI_Init(&argc,&argv);

  Homme::Context::singleton().create<Homme::Comm>().reset_mpi_comm(MPI_COMM_WORLD);
  Homme::initialize_hommexx_session();

  int result = Catch::Session().run(argc, argv);

  Homme::finalize_hommexx_session();

  // Finalize mpi
  MPI_Finalize();

  return result;
}
