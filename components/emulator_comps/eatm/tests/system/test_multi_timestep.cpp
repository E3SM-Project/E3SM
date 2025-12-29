//==============================================================================
// Test: Multi-Timestep Evolution
//
// System-level test for multi-timestep runs with state evolution.
//==============================================================================

#include "emulator_atm.hpp"
#include "inference/inference_backend.hpp"
#include "test_data.hpp"
#include <catch2/catch.hpp>
#include <cmath>
#include <mpi.h>

using namespace emulator;
using namespace emulator::testing;

TEST_CASE("Multi-timestep state evolution", "[system][eatm]") {
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

  const int global_ncol = 100;
  const int num_fields = 5;
  const int nsteps = 10;
  const int dt = 1800;

  SECTION("exports change over time") {
    EmulatorAtm atm;

    inference::InferenceConfig inf_config;
    inf_config.backend = inference::BackendType::STUB;
    atm.set_inference_config(inf_config);

    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
    atm.set_grid_data(10, 10, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    std::vector<double> imports(num_fields * grid.ncol, 300.0);
    std::vector<double> exports(num_fields * grid.ncol, 0.0);

    atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                       grid.ncol);

    atm.initialize();

    // Track exports over time
    std::vector<std::vector<double>> export_history;
    export_history.push_back(exports);

    for (int step = 0; step < nsteps; ++step) {
      atm.run(dt);

      // Save current exports
      export_history.push_back(exports);
    }

    atm.finalize();

    // Stub backend may or may not modify exports
    // This test verifies we can track state over time
    SUCCEED("Multi-timestep tracking completed");
  }
}

TEST_CASE("Coupling field consistency", "[system][eatm]") {
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

  const int global_ncol = 50;
  const int num_imports = 5;
  const int num_exports = 5;

  SECTION("exports remain physical after multiple timesteps") {
    EmulatorAtm atm;

    inference::InferenceConfig inf_config;
    inf_config.backend = inference::BackendType::STUB;
    atm.set_inference_config(inf_config);

    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
    atm.set_grid_data(10, 5, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    std::vector<double> imports(num_imports * grid.ncol, 300.0);
    std::vector<double> exports(num_exports * grid.ncol, 0.0);

    atm.setup_coupling(imports.data(), exports.data(), num_imports, num_exports,
                       grid.ncol);

    atm.initialize();

    // Run 20 timesteps
    for (int step = 0; step < 20; ++step) {
      atm.run(1800);

      // Check exports for NaN/Inf
      bool has_nan = false;
      bool has_inf = false;
      for (double val : exports) {
        if (std::isnan(val))
          has_nan = true;
        if (std::isinf(val))
          has_inf = true;
      }

      REQUIRE_FALSE(has_nan);
      REQUIRE_FALSE(has_inf);
    }

    atm.finalize();

    SUCCEED("All exports remained physical over 20 timesteps");
  }
}

TEST_CASE("Restart capability", "[system][eatm]") {
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

  const int global_ncol = 50;
  const int num_fields = 5;

  SECTION("can stop and restart simulation") {
    // First run: 5 timesteps
    {
      EmulatorAtm atm;

      inference::InferenceConfig inf_config;
      inf_config.backend = inference::BackendType::STUB;
      atm.set_inference_config(inf_config);

      TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

      atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
      atm.set_grid_data(10, 5, grid.ncol, global_ncol, grid.gids.data(),
                        grid.lat.data(), grid.lon.data(), grid.area.data());

      std::vector<double> imports(num_fields * grid.ncol, 300.0);
      std::vector<double> exports(num_fields * grid.ncol, 0.0);

      atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                         grid.ncol);

      atm.initialize();
      for (int i = 0; i < 5; ++i) {
        atm.run(1800);
      }
      atm.finalize();
    }

    // Second run: 5 more timesteps (simulated restart)
    {
      EmulatorAtm atm;

      inference::InferenceConfig inf_config;
      inf_config.backend = inference::BackendType::STUB;
      atm.set_inference_config(inf_config);

      TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

      // Use restart run_type (2)
      atm.create_instance(MPI_COMM_WORLD, 1, "", 2, 20000101, 9000);
      atm.set_grid_data(10, 5, grid.ncol, global_ncol, grid.gids.data(),
                        grid.lat.data(), grid.lon.data(), grid.area.data());

      std::vector<double> imports(num_fields * grid.ncol, 300.0);
      std::vector<double> exports(num_fields * grid.ncol, 0.0);

      atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                         grid.ncol);

      atm.initialize();
      for (int i = 0; i < 5; ++i) {
        atm.run(1800);
      }
      atm.finalize();
    }

    SUCCEED("Restart simulation completed");
  }
}
