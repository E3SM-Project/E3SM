// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("StubBackend factory creation", "[stub_backend]") {
  InferenceConfig config;
  config.input_channels = 4;
  config.output_channels = 2;

  auto backend = create_backend(BackendType::STUB, config);

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub");
}

TEST_CASE("StubBackend lifecycle", "[stub_backend]") {
  InferenceConfig config;
  config.input_channels = 4;
  config.output_channels = 2;

  auto backend = create_backend(BackendType::STUB, config);

  // Run inference (no-op: outputs unchanged)
  double inputs[4] = {1, 2, 3, 4};
  double outputs[2] = {99, 99};

  REQUIRE(backend->infer(inputs, outputs));

  REQUIRE(outputs[0] == 99.0);
  REQUIRE(outputs[1] == 99.0);

  // Finalize
  backend->finalize();
}

TEST_CASE("create_backend fallback for unknown type", "[stub_backend]") {
  InferenceConfig config;
  auto backend = create_backend(static_cast<BackendType>(999), config);

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub"); // Falls back to stub
}

} // namespace test
} // namespace inference
} // namespace emulator
