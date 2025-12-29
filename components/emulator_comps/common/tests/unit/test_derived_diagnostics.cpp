/**
 * @file test_derived_diagnostics.cpp
 * @brief Unit tests for derived diagnostics infrastructure.
 *
 * Tests pattern matching, factory creation, and computation.
 */

#include "diagnostics/diagnostic_factory.hpp"
#include "diagnostics/horiz_avg_diagnostic.hpp"
#include "diagnostics/vert_slice_diagnostic.hpp"
#include <catch2/catch.hpp>

using namespace emulator;
using Catch::Approx;

// =============================================================================
// Pattern Matching Tests
// =============================================================================

TEST_CASE("is_derived_diagnostic identifies patterns correctly",
          "[diagnostics][factory]") {

  SECTION("Horizontal average patterns") {
    CHECK(is_derived_diagnostic("temperature_horiz_avg"));
    CHECK(is_derived_diagnostic("surface_pressure_global_mean"));
    CHECK(is_derived_diagnostic("T_horiz_avg"));
  }

  SECTION("Vertical slice patterns") {
    CHECK(is_derived_diagnostic("air_temperature_at_lev0"));
    CHECK(is_derived_diagnostic("wind_at_lev5"));
    CHECK(is_derived_diagnostic("q_at_lev10"));
  }

  SECTION("Non-diagnostic patterns") {
    CHECK_FALSE(is_derived_diagnostic("temperature"));
    CHECK_FALSE(is_derived_diagnostic("surface_pressure"));
    CHECK_FALSE(is_derived_diagnostic("T"));
    CHECK_FALSE(is_derived_diagnostic("horiz_avg")); // No prefix
    CHECK_FALSE(is_derived_diagnostic("at_lev3"));   // No prefix
  }
}

TEST_CASE("get_base_field_name extracts correctly", "[diagnostics][factory]") {

  SECTION("Horizontal average") {
    CHECK(get_base_field_name("temperature_horiz_avg") == "temperature");
    CHECK(get_base_field_name("surface_pressure_global_mean") ==
          "surface_pressure");
  }

  SECTION("Vertical slice") {
    CHECK(get_base_field_name("air_temperature_at_lev0") == "air_temperature");
    CHECK(get_base_field_name("wind_at_lev5") == "wind");
  }

  SECTION("Non-diagnostic returns unchanged") {
    CHECK(get_base_field_name("temperature") == "temperature");
    CHECK(get_base_field_name("random_field") == "random_field");
  }
}

// =============================================================================
// Factory Tests
// =============================================================================

TEST_CASE("create_diagnostic returns correct types", "[diagnostics][factory]") {
  DiagnosticMetadata metadata;
  metadata.area_weights = {0.25, 0.25, 0.25, 0.25};
  metadata.nlevs = 8;

  SECTION("Creates HorizAvgDiagnostic") {
    auto diag = create_diagnostic("T_horiz_avg", metadata);
    REQUIRE(diag != nullptr);
    CHECK(diag->name() == "T_horiz_avg");
    CHECK(diag->source_field() == "T");
  }

  SECTION("Creates VertSliceDiagnostic") {
    auto diag = create_diagnostic("pressure_at_lev3", metadata);
    REQUIRE(diag != nullptr);
    CHECK(diag->name() == "pressure_at_lev3");
    CHECK(diag->source_field() == "pressure");
  }

  SECTION("Returns nullptr for non-diagnostic") {
    auto diag = create_diagnostic("temperature", metadata);
    CHECK(diag == nullptr);
  }
}

// =============================================================================
// Output Size Tests
// =============================================================================

TEST_CASE("Diagnostic output_size is correct", "[diagnostics]") {
  DiagnosticMetadata metadata;
  metadata.area_weights = {1.0, 1.0, 1.0, 1.0};
  metadata.nlevs = 5;

  SECTION("HorizAvgDiagnostic output size") {
    auto diag = create_diagnostic("field_horiz_avg", metadata);
    REQUIRE(diag != nullptr);
    // Returns nlevs (one average per level)
    CHECK(diag->output_size(100, 5) == 5);
    CHECK(diag->output_size(100, 1) == 1);
  }

  SECTION("VertSliceDiagnostic output size") {
    auto diag = create_diagnostic("field_at_lev2", metadata);
    REQUIRE(diag != nullptr);
    // Returns ncols (one value per column)
    CHECK(diag->output_size(100, 5) == 100);
    CHECK(diag->output_size(50, 10) == 50);
  }
}
