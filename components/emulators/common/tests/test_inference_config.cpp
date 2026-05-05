// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("InferenceConfig defaults", "[inference_config]") {
  InferenceConfig config;

  REQUIRE(config.input_channels == 0);
  REQUIRE(config.output_channels == 0);
  REQUIRE_FALSE(config.verbose);
}

TEST_CASE("InferenceConfig can be set", "[inference_config]") {
  InferenceConfig config;
  config.input_channels = 10;
  config.output_channels = 5;
  config.verbose = true;

  REQUIRE(config.input_channels == 10);
  REQUIRE(config.output_channels == 5);
  REQUIRE(config.verbose);
}

} // namespace test
} // namespace inference
} // namespace emulator
