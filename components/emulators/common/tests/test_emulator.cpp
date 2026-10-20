#include "test_emulator.hpp"
#include <emulator_c_api.hpp>
#include <emulator_config.hpp>
#include <inference_backend.hpp>
#include <stub_inference_backend.hpp>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <mpi.h>

namespace e3sm::emulator::tests {

TEST_CASE("Emulator initialization builds field layout") {

  const std::filesystem::path config_file =
      TEST_DATA_DIR "/simple_emulator_config.yaml";

  Emulator emu{
      CreateConfig{.path = std::move(config_file), .grid = make_test_grid()},
      std::unique_ptr<inference::InferenceBackend>(new inference::StubBackend(inference::InferenceConfig{})),
  };

  REQUIRE(emu.is_initialized());
  REQUIRE(emu.name() == "lnd-emulator");
  REQUIRE(emu.type() == emulator::EmulatorType::Component);

  REQUIRE(emu.field_descriptors().size() == 3);

  SECTION("temperature field") {
    const auto& field = emu.get_field_descriptor("temperature");

    REQUIRE(field.role == emulator::ModelFieldRole::Input);
    REQUIRE(field.size == 4);

    const std::vector<std::string> expected{"column"};
    REQUIRE(field.config.dimensions == expected);

    REQUIRE(emu.get_field("temperature").size() == 4);
  }

  SECTION("soil-moisture field") {
    const auto& field = emu.get_field_descriptor("soil_moisture");

    REQUIRE(field.role == emulator::ModelFieldRole::Input);
    REQUIRE(field.size == 12);

    const std::vector<std::string> expected{
        "column",
        "soil_levels",
    };
    REQUIRE(field.config.dimensions == expected);

    REQUIRE(emu.get_field("soil_moisture").size() == 12);
  }

  SECTION("surface-flux field") {
    const auto& field = emu.get_field_descriptor("surface_flux");

    REQUIRE(field.role == emulator::ModelFieldRole::Output);
    REQUIRE(field.size == 4);

    REQUIRE(emu.get_field("surface_flux").size() == 4);
  }

  SECTION("flat buffers have calculated sizes") {
    REQUIRE(emu.input_fields().size() == 16);
    REQUIRE(emu.output_fields().size() == 4);

    REQUIRE(std::ranges::all_of(emu.input_fields(),
                                [](double value) { return value == 0.0; }));

    REQUIRE(std::ranges::all_of(emu.output_fields(),
                                [](double value) { return value == 0.0; }));
  }
}

} // namespace e3sm::emulator::tests
