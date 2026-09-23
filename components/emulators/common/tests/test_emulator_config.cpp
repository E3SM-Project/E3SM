#define CATCH_CONFIG_MAIN
#include <algorithm>
#include <catch2/catch.hpp>
#include <emulator_config.hpp>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

TEST_CASE("Parse sample emulator config") {
  const std::filesystem::path config_file =
      TEST_DATA_DIR "/simple_emulator_config.yaml";

  const auto config = emulator::read_emulator_config(config_file);
  std::cout << config << std::endl;

  REQUIRE(config.name == "lnd-emulator");
  REQUIRE(config.type == "component");

  SECTION("model configuration") {
    REQUIRE(config.model.path == std::filesystem::path{"land_model.pt"});
    REQUIRE(config.model.input_layout == emulator::TensorLayout::SampleMajor);
    REQUIRE(config.model.output_layout == emulator::TensorLayout::FeatureMajor);
  }

  SECTION("sample-space dimensions") {
    const std::vector<std::string> expected_dimensions{
        "column",
        "soil_levels",
    };

    REQUIRE(config.dimensions == expected_dimensions);
  }

  SECTION("input fields") {
    REQUIRE(config.inputs.size() == 2);

    const auto temperature = std::ranges::find(config.inputs, "temperature",
                                               &emulator::FieldConfig::name);

    REQUIRE(temperature != config.inputs.end());

    const std::vector<std::string> expected_temperature_dimensions{
        "column",
    };

    REQUIRE(temperature->dimensions == expected_temperature_dimensions);

    const auto soil_moisture = std::ranges::find(config.inputs, "soil_moisture",
                                                 &emulator::FieldConfig::name);

    REQUIRE(soil_moisture != config.inputs.end());

    const std::vector<std::string> expected_soil_moisture_dimensions{
        "column",
        "soil_levels",
    };

    REQUIRE(soil_moisture->dimensions == expected_soil_moisture_dimensions);
  }

  SECTION("output fields") {
    REQUIRE(config.outputs.size() == 1);

    const auto surface_flux = std::ranges::find(config.outputs, "surface_flux",
                                                &emulator::FieldConfig::name);

    REQUIRE(surface_flux != config.outputs.end());

    const std::vector<std::string> expected_surface_flux_dimensions{
        "column",
    };

    REQUIRE(surface_flux->dimensions == expected_surface_flux_dimensions);
  }
}
