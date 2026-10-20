#include "inference_backend.hpp"
#include "inference_config.hpp"
#include "stub_inference_backend.hpp"
#include <memory>
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "test_emulator.hpp"
#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <filesystem>
#include <process.hpp>

using namespace e3sm::coupler;

namespace e3sm::emulator::tests {

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

  //FIXME: these are just dummy replacements for the former test emulator subclasses
  //FIXME: (see test_emulator.hpp)
  /*
  auto forcing_emulator = ForcingEmulator{create_forcing};
  auto test_emulator = Emulator{create_emu};
  */
  auto stub_forcing = std::unique_ptr<inference::InferenceBackend>(new inference::StubBackend(inference::InferenceConfig{}));
  auto forcing_emulator = Emulator{create_forcing, std::move(stub_forcing)};
  auto stub_test = std::unique_ptr<inference::InferenceBackend>(new inference::StubBackend(inference::InferenceConfig{}));
  auto test_emulator = Emulator{create_emu, std::move(stub_test)};
  REQUIRE(!forcing_emulator.input_fields().empty());
  REQUIRE(!test_emulator.input_fields().empty());
  REQUIRE(!forcing_emulator.output_fields().empty());
  REQUIRE(!test_emulator.output_fields().empty());

  // FIXME: see test_emulator.hpp
  /*
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
  */

  CouplerDriver driver;

  driver.add_component(forcing_emulator);
  driver.add_component(test_emulator);

  const std::string config = TEST_DATA_DIR "/simple_coupling_config.yaml";

  driver.initialize(config);
  const ProcessRunOpts opts{3600.0};
  driver.run_component(forcing_emulator.name(), opts);
  driver.run_component(test_emulator.name(), opts);

  /*
  REQUIRE(std::ranges::equal(forcing_emulator.temperature(),
                             test_emulator.temperature()));
  REQUIRE(std::ranges::equal(forcing_emulator.soil_moisture(),
                             test_emulator.soil_moisture()));

  REQUIRE_FALSE(std::ranges::equal(forcing_emulator.surface_flux(),
                                   test_emulator.surface_flux()));
  */
}

} // namespace e3sm::emulator::tests
