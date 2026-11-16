// moab_test_main.cpp
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <mpi.h>

#include <moab/Types.hpp>
#include <moab/iMOAB.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  const ErrCode init_error = iMOAB_Initialize(0, nullptr);
  if (init_error != moab::MB_SUCCESS) {
    MPI_Abort(MPI_COMM_WORLD, static_cast<int>(init_error));
    return 1;
  }

  int local_result = 0;
  {
    Catch::Session session;
    local_result = session.run(argc, argv);
  }

  const ErrCode finalize_error = iMOAB_Finalize();
  if (finalize_error != moab::MB_SUCCESS) {
    local_result = 1;
  }

  int global_result = 0;
  MPI_Allreduce(&local_result, &global_result, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);

  MPI_Finalize();
  return global_result;
}
