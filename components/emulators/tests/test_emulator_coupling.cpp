#ifndef EMULATOR_TEST_MPI
#define CATCH_CONFIG_MAIN
#endif
#include <catch2/catch.hpp>

#include "emulator_config.hpp"
#include "inference_backend.hpp"
#include "mpi.h"
#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <filesystem>
#include <memory>
#include <process.hpp>

namespace e3sm::emulator::tests {
using namespace coupler;
// An embedded interpreter does not see a virtual environment's packages,
//
// so the build passes the location of the ones it found
const std::string kPythonPath =
    std::string(EMULATOR_TEST_FIXTURE_DIR) + ":" + EMULATOR_TEST_SITE_PACKAGES;

static constexpr int GLOBAL_COLS = 4;
constexpr std::vector<int> gids(int size) {
  std::vector<int> result;
  result.reserve(size);
  for (int i = 0; i < size; ++i) {
    result[i] = i;
  }
  return result;
}

enum class TestModel {
  Forcing = 0,
  Land = 1,
};

inline GridDesc make_test_grid(MPI_Comm comm) {
  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  std::vector<int> local_col_extents;
  local_col_extents.resize(comm_size);
  // round robin decomp
  // for(int i=0;i<GLOBAL_COLS;)

  return GridDesc{
      .grid_type = 0,
      .nx = GLOBAL_COLS,
      .ny = 1,
      .num_global_columns = GLOBAL_COLS,
      .num_local_columns = -1,
      .comm = comm,
      .dimensions =
          {
              {
                  .name = "column",
                  .local_extent = 4,
                  .global_extent = GLOBAL_COLS,
              },
              {
                  .name = "soil_levels",
                  .local_extent = 3,
                  .global_extent = 3,
              },
          },
      .col_gids = gids(GLOBAL_COLS),
      .lat = {10.0, 20.0, 30.0, 40.0},
      .lon = {-100.0, -99.0, -98.0, -97.0},
      .area = {1.0, 1.0, 1.0, 1.0},
  };
}

inline EmulatorConfig make_test_config(std::string_view config_path,
                                       inference::BackendType type,
                                       TestModel model) {

  std::cout << "Reading " << config_path << std::endl;
  auto cfg = read_emulator_config(config_path);
  if (type == inference::BackendType::LIBTORCH) {
    return cfg;
  }
  cfg.model.options["python_path"] = kPythonPath;
  cfg.model.options["python_module"] =
      (model == TestModel::Land) ? "lnd_emulator" : "forcing_emulator";
  return cfg;
}

TEST_CASE("Emulator coupling across inference backends", "[coupling][mpi]") {

  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  REQUIRE(mpi_initialized != 0);

  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  MPI_Comm emu_comm = MPI_COMM_SELF;

  // MPI_Comm_split(emu_comm, (rank < 2) ? 0 : 1, rank, &emu_comm);

  const auto forcing_backend = GENERATE(EMULATOR_TEST_BACKENDS);
  const auto land_backend = GENERATE(EMULATOR_TEST_BACKENDS);

  CAPTURE(forcing_backend, land_backend);

  const auto grid = make_test_grid(emu_comm);

  auto forcing_config =
      make_test_config(TEST_DATA_DIR "/forcing_emulator_config.yaml",
                       forcing_backend, TestModel::Forcing);

  auto land_config =
      make_test_config(TEST_DATA_DIR "/land_emulator_config.yaml", land_backend,
                       TestModel::Land);

  std::cout << forcing_config << "\n" << land_config << std::endl;

  std::string forcing_name = forcing_config.name;
  std::string lnd_name = land_config.name;

  // Illustrative API: construct using the already-parsed configuration.

  CouplerDriver driver;
  driver.add_component(std::make_unique<Emulator>(grid, forcing_config));
  driver.add_component(std::make_unique<Emulator>(grid, land_config));
  driver.initialize(TEST_DATA_DIR "/test_coupling_plan.yaml");
  Emulator& forcing_emulator = driver.get_component_as<Emulator>(forcing_name);
  Emulator& land_emulator = driver.get_component_as<Emulator>(lnd_name);

  const ProcessRunOpts opts{3600.0};
  driver.run_component(forcing_emulator.name(), opts);
  driver.run_component(land_emulator.name(), opts);

  REQUIRE(std::ranges::equal(forcing_emulator.get_field("temperature"),
                             land_emulator.get_field("temperature")));

  REQUIRE(std::ranges::equal(forcing_emulator.get_field("soil_moisture"),
                             land_emulator.get_field("soil_moisture")));

  const auto temperature = land_emulator.get_field("temperature");
  const auto moisture = land_emulator.get_field("soil_moisture");
  const auto flux = land_emulator.get_field("surface_flux");

  REQUIRE(!temperature.empty());
  REQUIRE(moisture.size() % temperature.size() == 0);
  REQUIRE(flux.size() == temperature.size());

  const auto nlevels = moisture.size() / temperature.size();
  REQUIRE(nlevels > 0);

  for (std::size_t column = 0; column < temperature.size(); ++column) {
    double mean = 0.0;
    for (std::size_t level = 0; level < nlevels; ++level) {
      mean += moisture[column * nlevels + level];
    }
    mean /= static_cast<double>(nlevels);

    REQUIRE(flux[column] == Approx(2.0 * temperature[column] + 100.0 * mean));
  }
}
} // namespace e3sm::emulator::tests
