// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"
#include "stub_inference_backend.hpp"

#include <vector>

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("The factory returns a ready backend", "[stub_backend]") {
  InferenceConfig config;
  auto backend = create_backend(config, InferenceContext());

  REQUIRE(backend != nullptr);
  REQUIRE(backend->name() == "Stub");
  REQUIRE(backend->is_initialized());
}

TEST_CASE("The tensor path carries names and shapes", "[stub_backend]") {
  InferenceConfig config;
  config.backend = "stub";
  auto backend = create_backend(config, InferenceContext());

  std::vector<double> T(12, 300.0);
  std::vector<double> dT(12, -1.0);

  TensorMap inputs;
  inputs.wrap("T", static_cast<const double *>(T.data()), {4, 3});
  TensorMap outputs;
  outputs.wrap("dT", dT.data(), {4, 3});

  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(dT[0] == -1.0); // unchanged, so a forgotten field stays visible
}

TEST_CASE("The flat path wraps both buffers", "[stub_backend]") {
  InferenceConfig config;
  config.input_channels = 4;
  config.output_channels = 2;
  auto backend = create_backend(config, InferenceContext());

  const double in[4] = {1, 2, 3, 4};
  double out[2] = {99, 99};

  REQUIRE(backend->infer(in, out));
  REQUIRE(out[0] == 99.0); // the stub changes nothing
}

TEST_CASE("The flat path needs its channel counts", "[stub_backend]") {
  InferenceConfig config;
  auto backend = create_backend(config, InferenceContext());

  const double in[4] = {0, 0, 0, 0};
  double out[4] = {0, 0, 0, 0};
  REQUIRE_THROWS_AS(backend->infer(in, out, 1), InferenceError);
}

TEST_CASE("The lifecycle is idempotent at both ends", "[stub_backend]") {
  InferenceConfig config;
  auto backend = create_backend(config, InferenceContext());

  backend->initialize(); // already initialized by the factory
  REQUIRE(backend->is_initialized());

  // A component may finalize explicitly and then be destroyed.
  backend->finalize();
  backend->finalize();
  REQUIRE_FALSE(backend->is_initialized());
}

TEST_CASE("Inference before initialize() is refused", "[stub_backend]") {
  StubBackend backend{InferenceConfig(), InferenceContext()};

  TensorMap inputs;
  TensorMap outputs;
  REQUIRE_THROWS_AS(backend.infer(inputs, outputs), InferenceError);
}

TEST_CASE("An unknown backend name is refused", "[stub_backend]") {
  InferenceConfig config;
  config.backend = "tensorflow";
  REQUIRE_THROWS_WITH(create_backend(config, InferenceContext()),
                      Catch::Contains("This build has"));
}

TEST_CASE("available_backends always includes the stub", "[stub_backend]") {
  const auto names = available_backends();
  REQUIRE_FALSE(names.empty());
  REQUIRE(names.front() == "stub");
}

} // namespace test
} // namespace inference
} // namespace emulator
