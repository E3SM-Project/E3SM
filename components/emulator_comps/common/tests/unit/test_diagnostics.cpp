//==============================================================================
// Test: Diagnostics Configuration
//
// Unit tests for diagnostic output configuration and utilities.
//==============================================================================

#include "emulator_diagnostics.hpp"
#include <catch2/catch.hpp>

using namespace emulator;

// ============================================================================
// String Conversion Tests
// ============================================================================

TEST_CASE("FrequencyUnit string conversion", "[unit][diagnostics]") {
  SECTION("string to enum") {
    REQUIRE(str_to_freq_unit("nsteps") == FrequencyUnit::NSTEPS);
    REQUIRE(str_to_freq_unit("ndays") == FrequencyUnit::NDAYS);
    REQUIRE(str_to_freq_unit("NDAYS") == FrequencyUnit::NDAYS);
    REQUIRE(str_to_freq_unit("nmonths") == FrequencyUnit::NMONTHS);
    REQUIRE(str_to_freq_unit("nyears") == FrequencyUnit::NYEARS);
    REQUIRE(str_to_freq_unit("none") == FrequencyUnit::NONE);
    REQUIRE(str_to_freq_unit("invalid") == FrequencyUnit::NONE);
  }

  SECTION("enum to string") {
    REQUIRE(freq_unit_to_str(FrequencyUnit::NSTEPS) == "nsteps");
    REQUIRE(freq_unit_to_str(FrequencyUnit::NDAYS) == "ndays");
    REQUIRE(freq_unit_to_str(FrequencyUnit::NONE) == "none");
  }
}

TEST_CASE("OutputAvgType string conversion", "[unit][diagnostics]") {
  SECTION("string to enum") {
    REQUIRE(str_to_avg_type("instant") == OutputAvgType::INSTANT);
    REQUIRE(str_to_avg_type("INSTANT") == OutputAvgType::INSTANT);
    REQUIRE(str_to_avg_type("average") == OutputAvgType::AVERAGE);
    REQUIRE(str_to_avg_type("avg") == OutputAvgType::AVERAGE);
    REQUIRE(str_to_avg_type("mean") == OutputAvgType::AVERAGE);
    REQUIRE(str_to_avg_type("min") == OutputAvgType::MIN);
    REQUIRE(str_to_avg_type("max") == OutputAvgType::MAX);
    REQUIRE(str_to_avg_type("std") == OutputAvgType::STD);
    REQUIRE(str_to_avg_type("stddev") == OutputAvgType::STD);
    REQUIRE(str_to_avg_type("sum") == OutputAvgType::SUM);
    REQUIRE(str_to_avg_type("invalid") == OutputAvgType::INSTANT);
  }

  SECTION("enum to string") {
    REQUIRE(avg_type_to_str(OutputAvgType::INSTANT) == "instant");
    REQUIRE(avg_type_to_str(OutputAvgType::AVERAGE) == "average");
    REQUIRE(avg_type_to_str(OutputAvgType::STD) == "std");
    REQUIRE(avg_type_to_str(OutputAvgType::SUM) == "sum");
  }
}

TEST_CASE("OutputPrecision string conversion", "[unit][diagnostics]") {
  REQUIRE(str_to_precision("float32") == OutputPrecision::FLOAT32);
  REQUIRE(str_to_precision("float64") == OutputPrecision::FLOAT64);
  REQUIRE(str_to_precision("double") == OutputPrecision::FLOAT64);
  REQUIRE(str_to_precision("f64") == OutputPrecision::FLOAT64);
  REQUIRE(str_to_precision("invalid") == OutputPrecision::FLOAT32);

  REQUIRE(precision_to_str(OutputPrecision::FLOAT32) == "float32");
  REQUIRE(precision_to_str(OutputPrecision::FLOAT64) == "float64");
}

TEST_CASE("FileType suffix", "[unit][diagnostics]") {
  REQUIRE(file_type_suffix(FileType::HISTORY) == ".h.");
  REQUIRE(file_type_suffix(FileType::RESTART) == ".r.");
  REQUIRE(file_type_suffix(FileType::HISTORY_RESTART) == ".rh.");
}

// ============================================================================
// Configuration Structure Tests
// ============================================================================

TEST_CASE("OutputStreamConfig defaults", "[unit][diagnostics]") {
  OutputStreamConfig config;

  REQUIRE(config.stream_name == "h0");
  REQUIRE(config.filename_prefix == "emulator");
  REQUIRE(config.fields.empty());
  REQUIRE(config.frequency_unit == FrequencyUnit::NDAYS);
  REQUIRE(config.frequency == 1);
  REQUIRE(config.avg_type == OutputAvgType::INSTANT);
  REQUIRE(config.precision == OutputPrecision::FLOAT32);
  REQUIRE(config.max_snapshots_per_file == 1);
}

TEST_CASE("RestartConfig defaults", "[unit][diagnostics]") {
  RestartConfig config;

  REQUIRE(config.enabled == true);
  REQUIRE(config.filename_prefix == "emulator.r");
  REQUIRE(config.frequency_unit == FrequencyUnit::NDAYS);
  REQUIRE(config.frequency == 1);
}

TEST_CASE("DiagnosticConfig defaults", "[unit][diagnostics]") {
  DiagnosticConfig config;

  REQUIRE(config.history_streams.empty());
  REQUIRE(config.restart.enabled == true);
  REQUIRE(config.history_restart.enabled == true);
}
