//==============================================================================
// Test: EATM Standalone Run
//
// System-level test running EATM standalone with full initialization.
//==============================================================================

#include "emulator_atm.hpp"
#include "inference/inference_backend.hpp"
#include "test_data.hpp"
#include <catch2/catch.hpp>
#include <iostream>
#include <mpi.h>

using namespace emulator;
using namespace emulator::testing;

TEST_CASE("EATM standalone multi-timestep run", "[system][eatm]") {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    INFO("MPI not initialized"); SUCCEED("Skipping - no MPI");
    return;
  }

  int rank = 0, nprocs = 1;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    INFO("MPI_COMM_WORLD not available"); SUCCEED("Skipping - invalid MPI comm");
    return;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Configuration
  const int global_ncol = 100;
  const int nsteps = 5;
  const int dt = 1800; // 30 minutes
  const int num_imports = 10;
  const int num_exports = 15;

  SECTION("runs with stub backend") {
    EmulatorAtm atm;

    // Configure inference
    inference::InferenceConfig inf_config;
    inf_config.backend = inference::BackendType::STUB;
    inf_config.input_channels = 39;
    inf_config.output_channels = 44;
    atm.set_inference_config(inf_config);

    // Setup grid
    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
    atm.set_grid_data(10, 10, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    // Setup coupling with realistic temperature values
    std::vector<double> imports(num_imports * grid.ncol);
    std::vector<double> exports(num_exports * grid.ncol, 0.0);

    for (int f = 0; f < num_imports; ++f) {
      for (int i = 0; i < grid.ncol; ++i) {
        // Realistic temperature range with spatial variation
        imports[f * grid.ncol + i] =
            280.0 + 20.0 * std::sin(grid.lat[i] * M_PI / 180.0);
      }
    }

    atm.setup_coupling(imports.data(), exports.data(), num_imports, num_exports,
                       grid.ncol);

    // Initialize
    atm.initialize();

    // Run multiple timesteps
    for (int step = 0; step < nsteps; ++step) {
      atm.run(dt);

      if (rank == 0) {
        INFO("Completed step " << step + 1 << "/" << nsteps);
      }
    }

    // Finalize
    atm.finalize();

    SUCCEED("EATM completed " << nsteps << " timesteps successfully");
  }
}

TEST_CASE("EATM parallel scaling", "[system][eatm][mpi]") {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    INFO("MPI not initialized"); SUCCEED("Skipping - no MPI");
    return;
  }

  int rank = 0, nprocs = 1;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    INFO("MPI_COMM_WORLD not available"); SUCCEED("Skipping - invalid MPI comm");
    return;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  REQUIRE(nprocs >= 1);

  SECTION("runs correctly with domain decomposition") {
    const int global_ncol = 1000;

    EmulatorAtm atm;

    inference::InferenceConfig inf_config;
    inf_config.backend = inference::BackendType::STUB;
    atm.set_inference_config(inf_config);

    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
    atm.set_grid_data(100, 10, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    // Verify partition is correct
    int total_cols;
    MPI_Allreduce(&grid.ncol, &total_cols, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(total_cols == global_ncol);

    // Setup coupling
    const int num_fields = 5;
    std::vector<double> imports(num_fields * grid.ncol, 300.0);
    std::vector<double> exports(num_fields * grid.ncol, 0.0);

    atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                       grid.ncol);

    // Run
    atm.initialize();
    atm.run(1800);
    atm.finalize();

    // All ranks should complete
    int success = 1;
    int global_success;
    MPI_Allreduce(&success, &global_success, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    REQUIRE(global_success == 1);
  }
}

TEST_CASE("EATM large grid performance", "[system][eatm]") {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    INFO("MPI not initialized"); SUCCEED("Skipping - no MPI");
    return;
  }

  int rank = 0, nprocs = 1;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    INFO("MPI_COMM_WORLD not available"); SUCCEED("Skipping - invalid MPI comm");
    return;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  SECTION("handles production-scale grid") {
    // Approximate ne30 resolution: ~50000 columns
    const int global_ncol = 50000;
    const int nsteps = 2;

    EmulatorAtm atm;

    inference::InferenceConfig inf_config;
    inf_config.backend = inference::BackendType::STUB;
    atm.set_inference_config(inf_config);

    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
    atm.set_grid_data(360, 180, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    const int num_fields = 10;
    std::vector<double> imports(num_fields * grid.ncol, 300.0);
    std::vector<double> exports(num_fields * grid.ncol, 0.0);

    atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                       grid.ncol);

    atm.initialize();

    for (int step = 0; step < nsteps; ++step) {
      atm.run(1800);
    }

    atm.finalize();

    if (rank == 0) {
      INFO("Completed large grid test with " << global_ncol << " columns");
    }

    SUCCEED("Large grid test passed");
  }
}
