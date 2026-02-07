// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {
namespace test {

// ============================================================================
// BackendType Tests
// ============================================================================

TEST_CASE("backend_type_to_string", "[backend_type]") {
  REQUIRE(backend_type_to_string(BackendType::STUB) == "STUB");
}

TEST_CASE("backend_string_to_type", "[backend_type]") {
  REQUIRE(backend_string_to_type("STUB") == BackendType::STUB);
  REQUIRE(backend_string_to_type("stub") == BackendType::STUB);
  REQUIRE(backend_string_to_type("Stub") == BackendType::STUB);

  // Unknown falls back to STUB
  REQUIRE(backend_string_to_type("unknown") == BackendType::STUB);
  REQUIRE(backend_string_to_type("") == BackendType::STUB);
}

// ============================================================================
// InferenceConfig Tests
// ============================================================================

TEST_CASE("InferenceConfig defaults", "[inference_config]") {
  InferenceConfig config;

  REQUIRE(config.backend == BackendType::STUB);
  REQUIRE(config.model_path.empty());
  REQUIRE(config.device_id == -1);
  REQUIRE_FALSE(config.use_fp16);
  REQUIRE_FALSE(config.verbose);
  REQUIRE(config.input_channels == 44);
  REQUIRE(config.output_channels == 50);
  REQUIRE_FALSE(config.spatial_mode);
  REQUIRE_FALSE(config.dry_run);
}

// ============================================================================
// ValidationResult Tests
// ============================================================================

TEST_CASE("ValidationResult helpers", "[validation_result]") {
  ValidationResult result;

  REQUIRE(result.valid);
  REQUIRE(result.errors.empty());
  REQUIRE(result.warnings.empty());
  REQUIRE_FALSE(result.has_warnings());

  result.add_warning("This is a warning");
  REQUIRE(result.valid); // Warnings don't affect validity
  REQUIRE(result.has_warnings());
  REQUIRE(result.warnings.size() == 1);

  result.add_error("This is an error");
  REQUIRE_FALSE(result.valid);
  REQUIRE(result.errors.size() == 1);
}

} // namespace test
} // namespace inference
} // namespace emulator
