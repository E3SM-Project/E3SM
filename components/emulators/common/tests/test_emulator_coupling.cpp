#include <algorithm>
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "test_emulator.hpp"
#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <filesystem>
#include <process.hpp>

using namespace e3sm::coupler;

namespace emulator::tests {

static_assert(CoupledProcess<tests::ForcingEmulator>);
static_assert(CoupledProcess<tests::TestEmulator>);

TEST_CASE("Test Emulator Couplings") {

  const std::filesystem::path forcing_config_file =
      TEST_DATA_DIR "/forcing_emulator_config.yaml";
  const std::filesystem::path config_file =
      TEST_DATA_DIR "/simple_emulator_config.yaml";

  const auto grid = make_test_grid();
  CreateConfig create_forcing{.path = std::move(forcing_config_file),
                              .grid = std::move(grid)};

  CreateConfig create_emu{.path = std::move(config_file),
                          .grid = std::move(grid)};

  auto forcing_emulator = ForcingEmulator{create_forcing};
  auto test_emulator = TestEmulator{create_emu};
  forcing_emulator.init_data();
  test_emulator.init_data();
  REQUIRE(!forcing_emulator.input_fields().empty());
  REQUIRE(!test_emulator.input_fields().empty());
  REQUIRE(!forcing_emulator.output_fields().empty());
  REQUIRE(!test_emulator.output_fields().empty());

  const auto force_init_val = forcing_emulator.init_val;
  const auto test_init_val = test_emulator.init_val;

  REQUIRE(std::ranges::all_of(
      forcing_emulator.input_fields(),
      [force_init_val](double val) { return val == force_init_val; }));
  REQUIRE(std::ranges::all_of(
      test_emulator.input_fields(),
      [test_init_val](double val) { return val == test_init_val; }));

  REQUIRE(std::ranges::all_of(
      forcing_emulator.output_fields(),
      [force_init_val](double val) { return val == force_init_val; }));
  REQUIRE(std::ranges::all_of(
      test_emulator.output_fields(),
      [test_init_val](double val) { return val == test_init_val; }));

  CouplerDriver driver;

  driver.add_component(forcing_emulator);
  driver.add_component(test_emulator);

  const std::string config = TEST_DATA_DIR "/simple_coupling_config.yaml";

  driver.initialize(config);
  const ProcessRunOpts opts{3600.0};
  driver.run_component(forcing_emulator.name(), opts);
  driver.run_component(test_emulator.name(), opts);

  REQUIRE(std::ranges::equal(forcing_emulator.temperature(),
                             test_emulator.temperature()));
  REQUIRE(std::ranges::equal(forcing_emulator.soil_moisture(),
                             test_emulator.soil_moisture()));

  REQUIRE_FALSE(std::ranges::equal(forcing_emulator.surface_flux(),
                                   test_emulator.surface_flux()));
}

} // namespace emulator::tests
