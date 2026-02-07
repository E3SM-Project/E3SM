// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_factory.hpp"

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("StubBackend factory creation", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB);

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub");
  REQUIRE_FALSE(backend->is_initialized());
}

TEST_CASE("StubBackend lifecycle", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB);
  InferenceConfig config;
  config.input_channels = 4;
  config.output_channels = 2;

  // Initialize
  REQUIRE(backend->initialize(config));
  REQUIRE(backend->is_initialized());

  // Run inference (no-op: outputs unchanged)
  double inputs[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  double outputs[4] = {99, 99, 99, 99};

  REQUIRE(backend->infer(inputs, outputs));

  REQUIRE(outputs[0] == 99.0);
  REQUIRE(outputs[1] == 99.0);
  REQUIRE(outputs[2] == 99.0);
  REQUIRE(outputs[3] == 99.0);

  // Finalize
  backend->finalize();
  REQUIRE_FALSE(backend->is_initialized());
}

TEST_CASE("StubBackend infer fails when not initialized", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB);

  double inputs[4] = {1, 2, 3, 4};
  double outputs[2] = {0, 0};

  REQUIRE_FALSE(backend->infer(inputs, outputs));
}

TEST_CASE("create_backend fallback for unknown type", "[stub_backend]") {
  auto backend = create_backend(static_cast<BackendType>(999));

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub"); // Falls back to stub
}

} // namespace test
} // namespace inference
} // namespace emulator
