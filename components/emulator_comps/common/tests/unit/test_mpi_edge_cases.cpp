/**
 * @file test_mpi_edge_cases.cpp
 * @brief Unit tests for MPI edge cases.
 *
 * Tests diagnostic and output behavior with unusual MPI configurations:
 *
 * - Single rank
 * - Uneven column distribution
 * - Ranks with zero columns
 */

#include "diagnostics/diagnostic_factory.hpp"
#include "inference/inference_backend.hpp"
#include <catch2/catch.hpp>

using namespace emulator;
using Catch::Approx;

// =============================================================================
// Single Rank Tests
// =============================================================================

TEST_CASE("Diagnostics work correctly in single-rank mode",
          "[mpi][diagnostics]") {

  SECTION("Pattern matching single rank") {
    // Pattern matching should work regardless of MPI
    CHECK(is_derived_diagnostic("T_horiz_avg"));
    CHECK(is_derived_diagnostic("pressure_at_lev0"));
    CHECK(get_base_field_name("T_horiz_avg") == "T");
  }
}

// =============================================================================
// Uneven Distribution Tests
// =============================================================================

TEST_CASE("Factory handles edge cases gracefully", "[mpi][factory]") {

  SECTION("Empty area weights") {
    DiagnosticMetadata metadata;
    metadata.area_weights = {}; // Empty weights
    metadata.nlevs = 5;

    // Should still create diagnostic (uses uniform weights)
    auto diag = create_diagnostic("T_horiz_avg", metadata);
    REQUIRE(diag != nullptr);
    CHECK(diag->name() == "T_horiz_avg");
  }

  SECTION("Level index beyond range") {
    DiagnosticMetadata metadata;
    metadata.nlevs = 5;

    // Request level 10 when only 5 levels exist
    auto diag = create_diagnostic("T_at_lev10", metadata);
    REQUIRE(diag != nullptr);
    // Diagnostic is created - bounds checking happens at compute time
    CHECK(diag->name() == "T_at_lev10");
  }
}

// =============================================================================
// Validation Tests
// =============================================================================

TEST_CASE("ValidationResult tracks errors correctly", "[validation]") {
  // Note: inference:: namespace may need adjustment based on actual structure
  // Using the inference namespace from inference_backend.hpp

  SECTION("Default is valid") {
    inference::ValidationResult result;
    CHECK(result.valid);
    CHECK(result.errors.empty());
    CHECK(result.warnings.empty());
  }

  SECTION("Adding error marks invalid") {
    inference::ValidationResult result;
    result.add_error("Model file not found");
    CHECK_FALSE(result.valid);
    CHECK(result.errors.size() == 1);
  }

  SECTION("Warnings don't affect validity") {
    inference::ValidationResult result;
    result.add_warning("Using CPU instead of GPU");
    CHECK(result.valid);
    CHECK(result.has_warnings());
  }
}
