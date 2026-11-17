#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <mpi.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int local_result;
  {
    Catch::Session session;
    local_result = session.run(argc, argv);
  }

  int global_result = 0;
  MPI_Allreduce(
      &local_result, &global_result,
      1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  MPI_Finalize();
  return global_result;
}
