//==============================================================================
// Test: Inference Factory
//
// Unit tests for the inference backend factory.
//==============================================================================

#include "inference/inference_backend.hpp"
#include "inference/stub_backend.hpp"
#include <catch2/catch.hpp>

using namespace emulator::inference;

TEST_CASE("InferenceConfig defaults", "[unit][inference]") {
  InferenceConfig config;

  SECTION("has sensible defaults") {
    REQUIRE(config.backend == BackendType::STUB);
    REQUIRE(config.device_id == -1); // CPU
    REQUIRE(config.use_fp16 == false);
    REQUIRE(config.verbose == false);
    REQUIRE(config.input_channels == 44); // Default for ACE2-like models
    REQUIRE(config.output_channels == 50);
  }
}

TEST_CASE("Backend type conversion", "[unit][inference]") {
  SECTION("backend_type_to_string converts correctly") {
    REQUIRE(backend_type_to_string(BackendType::STUB) == "STUB");
    REQUIRE(backend_type_to_string(BackendType::LIBTORCH) == "LIBTORCH");
  }

  SECTION("parse_backend_type parses STUB variants") {
    REQUIRE(parse_backend_type("stub") == BackendType::STUB);
    REQUIRE(parse_backend_type("STUB") == BackendType::STUB);
    REQUIRE(parse_backend_type("unknown") == BackendType::STUB); // Default
  }

  SECTION("parse_backend_type parses LIBTORCH variants") {
    REQUIRE(parse_backend_type("LIBTORCH") == BackendType::LIBTORCH);
    REQUIRE(parse_backend_type("libtorch") == BackendType::LIBTORCH);
    REQUIRE(parse_backend_type("torch") == BackendType::LIBTORCH);
  }
}

TEST_CASE("create_backend factory function", "[unit][inference]") {
  SECTION("creates StubBackend for STUB type") {
    auto backend = create_backend(BackendType::STUB);

    REQUIRE(backend != nullptr);
    REQUIRE(backend->type() == BackendType::STUB);
    REQUIRE(backend->name() == "Stub");
  }

  SECTION("creates backend from config") {
    InferenceConfig config;
    config.backend = BackendType::STUB;
    config.input_channels = 20;
    config.output_channels = 10;

    auto backend = create_backend(config);

    REQUIRE(backend != nullptr);
    REQUIRE(backend->type() == BackendType::STUB);
  }

#ifdef EMULATOR_HAS_LIBTORCH
  SECTION("creates LibTorchBackend when available") {
    auto backend = create_backend(BackendType::LIBTORCH);
    REQUIRE(backend != nullptr);
    REQUIRE(backend->type() == BackendType::LIBTORCH);
  }
#endif
}

TEST_CASE("Backend lifecycle", "[unit][inference]") {
  auto backend = create_backend(BackendType::STUB);
  REQUIRE(backend != nullptr);

  SECTION("full lifecycle: initialize -> infer -> finalize") {
    InferenceConfig config;
    config.input_channels = 5;
    config.output_channels = 3;

    // Initialize
    REQUIRE(backend->initialize(config));
    REQUIRE(backend->is_initialized());

    // Infer
    std::vector<double> inputs(5, 1.0);
    std::vector<double> outputs(3, 0.0);
    REQUIRE(backend->infer(inputs.data(), outputs.data(), 1));

    // Finalize
    backend->finalize();
    REQUIRE_FALSE(backend->is_initialized());
  }
}
