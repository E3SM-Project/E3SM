//==============================================================================
// Test: EmulatorAtm Component
//
// Component-level tests for the atmosphere emulator.
//==============================================================================

#include "emulator_atm.hpp"
#include "inference/inference_backend.hpp"
#include "test_data.hpp"
#include <catch2/catch.hpp>
#include <mpi.h>

using namespace emulator;
using namespace emulator::testing;

TEST_CASE("EmulatorAtm construction", "[component][eatm]") {
  EmulatorAtm atm;

  SECTION("can be default constructed") {
    REQUIRE(atm.type() == CompType::ATM);
  }
}

TEST_CASE("EmulatorAtm inference configuration", "[component][eatm]") {
  EmulatorAtm atm;

  SECTION("can set and get inference config") {
    inference::InferenceConfig config;
    config.backend = inference::BackendType::STUB;
    config.input_channels = 39;
    config.output_channels = 44;

    atm.set_inference_config(config);

    auto retrieved = atm.get_inference_config();
    REQUIRE(retrieved.backend == inference::BackendType::STUB);
    REQUIRE(retrieved.input_channels == 39);
    REQUIRE(retrieved.output_channels == 44);
  }
}

TEST_CASE("EmulatorAtm initialization with grid", "[component][eatm]") {
  // Check if MPI is functional
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

  EmulatorAtm atm;

  // Configure inference backend
  inference::InferenceConfig inf_config;
  inf_config.backend = inference::BackendType::STUB;
  atm.set_inference_config(inf_config);

  SECTION("can be initialized with grid data") {
    const int global_ncol = 100;
    TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

    // Create instance
    atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);

    // Set grid data
    atm.set_grid_data(10, 10, grid.ncol, global_ncol, grid.gids.data(),
                      grid.lat.data(), grid.lon.data(), grid.area.data());

    REQUIRE(atm.get_num_local_cols() == grid.ncol);
    REQUIRE(atm.get_num_global_cols() == global_ncol);
  }
}

TEST_CASE("EmulatorAtm coupling setup", "[component][eatm]") {
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

  EmulatorAtm atm;

  // Configure
  inference::InferenceConfig inf_config;
  inf_config.backend = inference::BackendType::STUB;
  atm.set_inference_config(inf_config);

  const int global_ncol = 50;
  TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

  atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
  atm.set_grid_data(10, 5, grid.ncol, global_ncol, grid.gids.data(),
                    grid.lat.data(), grid.lon.data(), grid.area.data());

  SECTION("can setup coupling buffers") {
    const int num_imports = 10;
    const int num_exports = 15;

    std::vector<double> imports(num_imports * grid.ncol, 0.0);
    std::vector<double> exports(num_exports * grid.ncol, 0.0);

    // Fill imports with test data
    for (int f = 0; f < num_imports; ++f) {
      for (int i = 0; i < grid.ncol; ++i) {
        imports[f * grid.ncol + i] = 300.0 + f * 10.0; // Temperature-like
      }
    }

    atm.setup_coupling(imports.data(), exports.data(), num_imports, num_exports,
                       grid.ncol);

    // If we get here without crash, coupling setup worked
    SUCCEED("Coupling buffers setup successfully");
  }
}

TEST_CASE("EmulatorAtm single timestep", "[component][eatm]") {
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

  EmulatorAtm atm;

  // Configure
  inference::InferenceConfig inf_config;
  inf_config.backend = inference::BackendType::STUB;
  atm.set_inference_config(inf_config);

  const int global_ncol = 100;
  const int num_imports = 5;
  const int num_exports = 5;

  TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

  atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
  atm.set_grid_data(10, 10, grid.ncol, global_ncol, grid.gids.data(),
                    grid.lat.data(), grid.lon.data(), grid.area.data());

  std::vector<double> imports(num_imports * grid.ncol, 300.0);
  std::vector<double> exports(num_exports * grid.ncol, 0.0);

  atm.setup_coupling(imports.data(), exports.data(), num_imports, num_exports,
                     grid.ncol);

  SECTION("can run initialize-run-finalize cycle") {
    atm.initialize();
    atm.run(1800); // 30 minute timestep
    atm.finalize();

    SUCCEED("Single timestep completed successfully");
  }
}

TEST_CASE("EmulatorAtm parallel consistency", "[component][eatm][mpi]") {
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

  EmulatorAtm atm;

  inference::InferenceConfig inf_config;
  inf_config.backend = inference::BackendType::STUB;
  atm.set_inference_config(inf_config);

  const int global_ncol = 200;
  TestGrid grid = create_partitioned_grid(global_ncol, rank, nprocs);

  atm.create_instance(MPI_COMM_WORLD, 1, "", 1, 20000101, 0);
  atm.set_grid_data(20, 10, grid.ncol, global_ncol, grid.gids.data(),
                    grid.lat.data(), grid.lon.data(), grid.area.data());

  SECTION("partition covers all global columns") {
    int total_cols;
    MPI_Allreduce(&grid.ncol, &total_cols, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(total_cols == global_ncol);
  }

  SECTION("all ranks can complete run") {
    const int num_fields = 5;
    std::vector<double> imports(num_fields * grid.ncol, 300.0);
    std::vector<double> exports(num_fields * grid.ncol, 0.0);

    atm.setup_coupling(imports.data(), exports.data(), num_fields, num_fields,
                       grid.ncol);

    atm.initialize();
    atm.run(1800);
    atm.finalize();

    // Synchronize and check all ranks succeeded
    int local_success = 1;
    int global_success;
    MPI_Allreduce(&local_success, &global_success, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    REQUIRE(global_success == 1);
  }
}
