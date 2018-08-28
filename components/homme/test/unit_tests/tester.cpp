
#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include <Kokkos_Core.hpp>

#include <Hommexx_Session.hpp>

#ifndef HOMMEXX_NO_MPI_CONTEXT
  #include <mpi/Comm.hpp>
  #include <mpi/MpiContext.hpp>
#endif

#include <mpi.h>

int main(int argc, char **argv) {

  // Initialize mpi
  MPI_Init(&argc,&argv);

#ifndef HOMMEXX_NO_MPI_CONTEXT
  Homme::MpiContext::singleton().get_comm().reset_mpi_comm(MPI_COMM_WORLD);
#endif
  Homme::initialize_hommexx_session();

  int result = Catch::Session().run(argc, argv);

  Homme::finalize_hommexx_session();

  // Finalize mpi
  MPI_Finalize();

  return result;
}
