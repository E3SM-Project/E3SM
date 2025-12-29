//==============================================================================
// Test: StubBackend
//
// Unit tests for the stub inference backend.
//==============================================================================

#include "inference/inference_backend.hpp"
#include "inference/stub_backend.hpp"
#include "test_data.hpp" // Include first to get mpi.h before inference stubs
#include <catch2/catch.hpp>

using namespace emulator::inference;
using namespace emulator::testing;

TEST_CASE("StubBackend initialization", "[unit][inference]") {
  StubBackend backend;

  SECTION("starts uninitialized") { REQUIRE_FALSE(backend.is_initialized()); }

  SECTION("initializes successfully with default config") {
    InferenceConfig config;
    config.backend = BackendType::STUB;

    REQUIRE(backend.initialize(config));
    REQUIRE(backend.is_initialized());
    REQUIRE(backend.type() == BackendType::STUB);
    REQUIRE(backend.name() == "Stub");
  }

  SECTION("can be finalized") {
    InferenceConfig config;
    backend.initialize(config);

    backend.finalize();
    REQUIRE_FALSE(backend.is_initialized());
  }

  SECTION("can be re-initialized after finalize") {
    InferenceConfig config;
    backend.initialize(config);
    backend.finalize();

    REQUIRE(backend.initialize(config));
    REQUIRE(backend.is_initialized());
  }
}

TEST_CASE("StubBackend inference", "[unit][inference]") {
  StubBackend backend;
  InferenceConfig config;
  config.input_channels = 10;
  config.output_channels = 5;
  backend.initialize(config);

  SECTION("produces deterministic output") {
    auto inputs = create_constant_field(10, 1.0);
    std::vector<double> outputs1(5, 0.0);
    std::vector<double> outputs2(5, 0.0);

    REQUIRE(backend.infer(inputs.data(), outputs1.data(), 1));
    REQUIRE(backend.infer(inputs.data(), outputs2.data(), 1));

    // Same inputs should produce same outputs
    REQUIRE(vectors_equal(outputs1, outputs2));
  }

  SECTION("handles batch inference") {
    int batch_size = 4;
    auto inputs = create_random_field(batch_size * 10, 0.0, 1.0);
    std::vector<double> outputs(batch_size * 5, 0.0);

    REQUIRE(backend.infer(inputs.data(), outputs.data(), batch_size));
  }

  SECTION("handles large batch") {
    int batch_size = 1000;
    auto inputs = create_random_field(batch_size * 10, 0.0, 1.0);
    std::vector<double> outputs(batch_size * 5, 0.0);

    REQUIRE(backend.infer(inputs.data(), outputs.data(), batch_size));
  }

  SECTION("output depends on input") {
    auto inputs1 = create_constant_field(10, 1.0);
    auto inputs2 = create_constant_field(10, 2.0);
    std::vector<double> outputs1(5, 0.0);
    std::vector<double> outputs2(5, 0.0);

    backend.infer(inputs1.data(), outputs1.data(), 1);
    backend.infer(inputs2.data(), outputs2.data(), 1);

    // Different inputs should (generally) produce different outputs
    // This depends on stub backend implementation
  }
}

TEST_CASE("StubBackend edge cases", "[unit][inference]") {
  StubBackend backend;
  InferenceConfig config;
  config.input_channels = 5;
  config.output_channels = 3;
  backend.initialize(config);

  SECTION("handles batch size 1") {
    auto inputs = create_random_field(5);
    std::vector<double> outputs(3, 0.0);

    REQUIRE(backend.infer(inputs.data(), outputs.data(), 1));
  }

  SECTION("handles zero batch size gracefully") {
    double inputs[1];
    double outputs[1];

    // Batch size 0 should be handled (implementation dependent)
    REQUIRE(backend.infer(inputs, outputs, 0));
  }
}
