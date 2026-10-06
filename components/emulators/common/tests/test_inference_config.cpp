// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("InferenceConfig defaults", "[inference_config]") {
  InferenceConfig config;

  REQUIRE(config.model_path.empty());
  REQUIRE(config.input_channels == 0);
  REQUIRE(config.output_channels == 0);
  REQUIRE_FALSE(config.verbose);
  REQUIRE(config.options.empty());
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

TEST_CASE("InferenceConfig options fall back when absent",
          "[inference_config]") {
  InferenceConfig config;

  REQUIRE(config.get("device") == "");
  REQUIRE(config.get("device", "cpu") == "cpu");
  REQUIRE(config.get_int("num_threads", 4) == 4);
  REQUIRE(config.get_bool("jit_optimize", true));
}

TEST_CASE("InferenceConfig options are typed on read", "[inference_config]") {
  InferenceConfig config;
  config.set("device", "cuda:1");
  config.set("num_threads", "8");
  config.set("jit_optimize", "true");

  REQUIRE(config.get("device", "cpu") == "cuda:1");
  REQUIRE(config.get_int("num_threads", 0) == 8);
  REQUIRE(config.get_bool("jit_optimize", false));

  config.set("jit_optimize", "0");
  REQUIRE_FALSE(config.get_bool("jit_optimize", true));
}

TEST_CASE("InferenceConfig rejects malformed options", "[inference_config]") {
  InferenceConfig config;
  config.set("num_threads", "8x");
  config.set("jit_optimize", "maybe");

  REQUIRE_THROWS_AS(config.get_int("num_threads", 0), InferenceError);
  REQUIRE_THROWS_AS(config.get_bool("jit_optimize", false), InferenceError);

  config.set("num_threads", "");
  REQUIRE_THROWS_AS(config.get_int("num_threads", 0), InferenceError);
}

} // namespace test
} // namespace inference
} // namespace emulator
