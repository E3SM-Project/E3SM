// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"

#include <vector>

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

TEST_CASE("StubBackend tensor inference", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB, InferenceConfig());

  const std::vector<double> T(12, 300.0);
  std::vector<double> dT(12, -1.0);

  TensorMap inputs;
  inputs.wrap("T", T.data(), {4, 3});
  TensorMap outputs;
  outputs.wrap("dT", dT.data(), {4, 3});

  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(dT[0] == -1.0);
}

TEST_CASE("Flat-array inference needs channel counts", "[stub_backend]") {
  auto backend = create_backend(BackendType::STUB, InferenceConfig());

  const double inputs[4] = {1, 2, 3, 4};
  double outputs[2] = {99, 99};

  REQUIRE_THROWS_AS(backend->infer(inputs, outputs), InferenceError);
}

TEST_CASE("create_backend refuses an unknown type", "[stub_backend]") {
  InferenceConfig config;

  REQUIRE_THROWS_AS(create_backend(static_cast<BackendType>(999), config),
                    InferenceError);
}

} // namespace test
} // namespace inference
} // namespace emulator
