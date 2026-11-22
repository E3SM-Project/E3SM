#ifndef EMULATOR_TEST_MPI
#define CATCH_CONFIG_MAIN
#include "field_types.hpp"
#endif
#include <catch2/catch.hpp>

#include "emulator_config.hpp"
#include "inference_backend.hpp"
#include "mpi.h"
#include <algorithm>
#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <filesystem>
#include <memory>
#include <process.hpp>
#include <ranges>

namespace e3sm::emulator::tests {
using namespace coupler;
using Integer = std::int32_t;
static constexpr double MASK_VAL = -9.99;

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

// Returns expected result vector that is MASK_VAL for inactive columns
std::vector<double> flux_expected(std::span<Integer> filter,
                                  ConstFieldView temperature,
                                  ConstFieldView moisture) {
  const auto nlevels = moisture.size() / temperature.size();
  REQUIRE(nlevels > 0);
  auto is_active = [](Integer c) { return c != -1; };
  auto flux_fn = [&](Integer c, double mean) {
    return 2.0 * temperature.as<double>()[c] + 100.0 * mean;
  };
  std::size_t num_active = std::ranges::count_if(filter, is_active);
  std::vector<double> result;
  result.resize(temperature.size(), MASK_VAL);

  for (Integer fc : filter | std::views::filter(is_active)) {
    double mean = 0.0;
    for (std::size_t level = 0; level < nlevels; ++level) {
      mean += moisture.as<double>()[fc * nlevels + level];
    }
    mean /= static_cast<double>(nlevels);
    result[fc] = flux_fn(fc, mean);
  }

  REQUIRE(std::ranges::count_if(result, [](double val) {
            return val != MASK_VAL;
          }) == num_active);

  return result;
}

TEST_CASE("Emulator coupling across inference backends", "[coupling][mpi]") {

  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  REQUIRE(mpi_initialized != 0);

  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  MPI_Comm emu_comm = MPI_COMM_SELF;

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

  auto coupled_fields = coupler::read_coupling_fields_from_yaml(
      TEST_DATA_DIR "/test_coupling_plan.yaml");

  CouplerDriver driver;
  driver.add_component(std::make_unique<Emulator>(grid, forcing_config));
  driver.add_component(std::make_unique<Emulator>(grid, land_config));
  driver.initialize(coupled_fields);
  Emulator& forcing_emulator = driver.get_component_as<Emulator>(forcing_name);
  Emulator& land_emulator = driver.get_component_as<Emulator>(lnd_name);

  ConstFieldView temp_forcing =
      forcing_emulator.get_field("temperature", ModelFieldRole::Output);
  ConstFieldView moisture_forcing =
      forcing_emulator.get_field("soil_moisture", ModelFieldRole::Output);

  const auto temperature_lnd =
      land_emulator.get_field("temperature", ModelFieldRole::Input);
  const auto moisture_lnd =
      land_emulator.get_field("soil_moisture", ModelFieldRole::Input);
  const auto flux =
      land_emulator.get_field("surface_flux", ModelFieldRole::Output);

  std::span<Integer> filter_column =
      land_emulator.get_field("filter_columns", ModelFieldRole::Input)
          .as<Integer>();

  constexpr auto cols = [] {
    std::array<Integer, GLOBAL_COLS> values{};
    std::iota(values.begin(), values.end(), 0);
    return values;
  }();

  const ProcessRunOpts opts{3600.0};
  // Filter columns in this case is NOT a coupled field. There we must supply it
  // explicitly in this test with each time step
  //
  for (int step = 0; step < GLOBAL_COLS; ++step) {

    // Arbitrary manner to change the filter each step. Example for cols =
    // {0,1,2,3}
    //   step N = filter
    //        0 = {0,-1,-1,-1}
    //        1 = {0, -1, -1, 3}
    //        2 = {0, -1, 2, -1}
    //        3 = {0, 1, 2, 3}
    std::ranges::transform(cols, filter_column.begin(), [step](Integer c) {
      return (c % (GLOBAL_COLS - step) == 0) ? c : -1;
    });

    driver.run_component(forcing_emulator.name(), opts);
    driver.run_component(land_emulator.name(), opts);

    REQUIRE(std::ranges::equal(temp_forcing.as<double>(),
                               temperature_lnd.as<double>()));

    REQUIRE(std::ranges::equal(moisture_forcing.as<double>(),
                               moisture_forcing.as<double>()));

    REQUIRE(!temperature_lnd.empty());
    REQUIRE(moisture_lnd.size() % temperature_lnd.size() == 0);
    REQUIRE(flux.size() == temperature_lnd.size());

    std::vector<double> expected_flux =
        flux_expected(filter_column, temperature_lnd, moisture_lnd);
    for (Integer col : cols) {
      REQUIRE(flux.as<double>()[col] == expected_flux[col]);
    }
  }
}

} // namespace e3sm::emulator::tests
