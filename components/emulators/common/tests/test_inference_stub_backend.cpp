// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_backend.hpp"
#include "stub_backend.hpp"

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("StubBackend factory creation", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB);

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub");
  REQUIRE(backend->type() == BackendType::STUB);
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

  // Run inference
  double inputs[8] = {1, 2, 3, 4, 5, 6, 7, 8}; // 2 grid points * 4 channels
  double outputs[4] = {99, 99, 99, 99};        // 2 grid points * 2 channels

  REQUIRE(backend->infer(inputs, outputs, 2));

  // Stub zeros outputs
  REQUIRE(outputs[0] == 0.0);
  REQUIRE(outputs[1] == 0.0);
  REQUIRE(outputs[2] == 0.0);
  REQUIRE(outputs[3] == 0.0);

  // Finalize
  backend->finalize();
  REQUIRE_FALSE(backend->is_initialized());
}

TEST_CASE("StubBackend infer fails when not initialized", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB);

  double inputs[4] = {1, 2, 3, 4};
  double outputs[2] = {0, 0};

  REQUIRE_FALSE(backend->infer(inputs, outputs, 1));
}

TEST_CASE("create_backend fallback for unknown type", "[stub_backend]") {
  // Use an invalid enum value cast (for testing purposes)
  auto backend = create_backend(static_cast<BackendType>(999));

  REQUIRE(backend != nullptr);
  REQUIRE(backend->type() == BackendType::STUB); // Falls back to stub
}

} // namespace test
} // namespace inference
} // namespace emulator
