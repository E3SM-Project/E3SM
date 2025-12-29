//==============================================================================
// Test: EmulatorConfig
//
// Unit tests for YAML configuration parsing.
//==============================================================================

#include "emulator_config.hpp"
#include "test_config.hpp"
#include <catch2/catch.hpp>

#include <cstdio>
#include <fstream>

using namespace emulator;

TEST_CASE("EmulatorConfig default values", "[unit][config]") {
  EmulatorConfig config;

  SECTION("has sensible defaults") {
    REQUIRE(config.runtime.enabled == true);
    REQUIRE(config.build.inference_backend == "stub");
    REQUIRE(config.coupling.zero_init_exports == true);
    REQUIRE(config.coupling.debug == false);
  }
}

TEST_CASE("EmulatorConfig YAML parsing", "[unit][config]") {
  // Create a temporary YAML file for testing
  const char *test_yaml = "test_config_temp.yaml";

  SECTION("parses valid YAML configuration") {
    // Write test YAML with new structure
    std::ofstream ofs(test_yaml);
    ofs << "eatm:\n"
        << "  build:\n"
        << "    grid_name: \"gauss180x360\"\n"
        << "    inference_backend: \"stub\"\n"
        << "  runtime:\n"
        << "    model_path: \"/path/to/model.pt\"\n"
        << "    ic_file: \"/path/to/ic.nc\"\n"
        << "    enabled: true\n"
        << "  coupling:\n"
        << "    zero_init_exports: true\n"
        << "    debug: false\n";
    ofs.close();

    auto config = parse_emulator_config(test_yaml, "eatm");

    REQUIRE(config.build.grid_name == "gauss180x360");
    REQUIRE(config.build.inference_backend == "stub");
    REQUIRE(config.runtime.model_path == "/path/to/model.pt");
    REQUIRE(config.runtime.ic_file == "/path/to/ic.nc");
    REQUIRE(config.runtime.enabled == true);
    REQUIRE(config.coupling.zero_init_exports == true);
    REQUIRE(config.coupling.debug == false);

    std::remove(test_yaml);
  }

  SECTION("parses libtorch backend configuration") {
    std::ofstream ofs(test_yaml);
    ofs << "eatm:\n"
        << "  build:\n"
        << "    inference_backend: \"libtorch\"\n"
        << "  runtime:\n"
        << "    model_path: \"/models/ace.pt\"\n"
        << "    ic_file: \"/data/ic.nc\"\n";
    ofs.close();

    auto config = parse_emulator_config(test_yaml, "eatm");

    REQUIRE(config.build.inference_backend == "libtorch");
    REQUIRE(config.runtime.model_path == "/models/ace.pt");
    REQUIRE(config.runtime.ic_file == "/data/ic.nc");

    std::remove(test_yaml);
  }

  SECTION("parses debug mode") {
    std::ofstream ofs(test_yaml);
    ofs << "eatm:\n"
        << "  coupling:\n"
        << "    debug: true\n";
    ofs.close();

    auto config = parse_emulator_config(test_yaml, "eatm");

    REQUIRE(config.coupling.debug == true);

    std::remove(test_yaml);
  }
}

TEST_CASE("EmulatorConfig with defaults fallback", "[unit][config]") {
  SECTION("returns defaults for nonexistent file") {
    auto config = parse_emulator_config_with_defaults(
        "nonexistent_file_12345.yaml", "eatm", false);

    REQUIRE(config.runtime.enabled == true);
    REQUIRE(config.build.inference_backend == "stub");
  }

  SECTION("returns defaults for missing section") {
    const char *test_yaml = "test_config_temp2.yaml";

    std::ofstream ofs(test_yaml);
    ofs << "other_section:\n"
        << "  some_key: some_value\n";
    ofs.close();

    auto config = parse_emulator_config_with_defaults(test_yaml, "eatm", false);

    // Should get defaults since 'eatm' section doesn't exist
    REQUIRE(config.runtime.enabled == true);
    REQUIRE(config.build.inference_backend == "stub");

    std::remove(test_yaml);
  }
}
