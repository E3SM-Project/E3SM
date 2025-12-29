//==============================================================================
// Test: Output Averaging
//
// Unit tests for output stream averaging logic.
//==============================================================================

#include "emulator_diagnostics.hpp"
#include <algorithm>
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

using namespace emulator;

// Helper function that mimics the averaging logic in EmulatorOutputStream
std::vector<double>
compute_average(const std::vector<std::vector<double>> &samples,
                OutputAvgType avg_type) {
  if (samples.empty()) {
    return {};
  }

  const int nsamples = static_cast<int>(samples.size());
  const int size = static_cast<int>(samples[0].size());
  std::vector<double> buffer(size, 0.0);
  std::vector<double> buffer_sq(size, 0.0);

  for (int s = 0; s < nsamples; ++s) {
    const auto &sample = samples[s];

    switch (avg_type) {
    case OutputAvgType::AVERAGE:
    case OutputAvgType::SUM:
      for (int i = 0; i < size; ++i) {
        buffer[i] += sample[i];
      }
      break;

    case OutputAvgType::MIN:
      if (s == 0) {
        buffer = sample;
      } else {
        for (int i = 0; i < size; ++i) {
          buffer[i] = std::min(buffer[i], sample[i]);
        }
      }
      break;

    case OutputAvgType::MAX:
      if (s == 0) {
        buffer = sample;
      } else {
        for (int i = 0; i < size; ++i) {
          buffer[i] = std::max(buffer[i], sample[i]);
        }
      }
      break;

    case OutputAvgType::STD:
      for (int i = 0; i < size; ++i) {
        buffer[i] += sample[i];
        buffer_sq[i] += sample[i] * sample[i];
      }
      break;

    default:
      break;
    }
  }

  // Compute final values
  std::vector<double> output(size);

  switch (avg_type) {
  case OutputAvgType::AVERAGE:
    for (int i = 0; i < size; ++i) {
      output[i] = buffer[i] / nsamples;
    }
    break;

  case OutputAvgType::SUM:
  case OutputAvgType::MIN:
  case OutputAvgType::MAX:
    output = buffer;
    break;

  case OutputAvgType::STD:
    for (int i = 0; i < size; ++i) {
      double mean = buffer[i] / nsamples;
      double mean_sq = buffer_sq[i] / nsamples;
      double variance = mean_sq - mean * mean;
      output[i] = std::sqrt(std::max(0.0, variance));
    }
    break;

  default:
    output = buffer;
    break;
  }

  return output;
}

TEST_CASE("Average single sample", "[unit][diagnostics]") {
  std::vector<std::vector<double>> samples = {{1.0, 2.0, 3.0}};

  auto result = compute_average(samples, OutputAvgType::AVERAGE);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(1.0));
  REQUIRE(result[1] == Approx(2.0));
  REQUIRE(result[2] == Approx(3.0));
}

TEST_CASE("Average multiple samples", "[unit][diagnostics]") {
  std::vector<std::vector<double>> samples = {
      {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {2.0, 3.0, 4.0}};

  auto result = compute_average(samples, OutputAvgType::AVERAGE);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(2.0));
  REQUIRE(result[1] == Approx(3.0));
  REQUIRE(result[2] == Approx(4.0));
}

TEST_CASE("Min multiple samples", "[unit][diagnostics]") {
  std::vector<std::vector<double>> samples = {
      {3.0, 4.0, 5.0}, {1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}};

  auto result = compute_average(samples, OutputAvgType::MIN);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(1.0));
  REQUIRE(result[1] == Approx(2.0));
  REQUIRE(result[2] == Approx(3.0));
}

TEST_CASE("Max multiple samples", "[unit][diagnostics]") {
  std::vector<std::vector<double>> samples = {
      {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {2.0, 3.0, 4.0}};

  auto result = compute_average(samples, OutputAvgType::MAX);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(3.0));
  REQUIRE(result[1] == Approx(4.0));
  REQUIRE(result[2] == Approx(5.0));
}

TEST_CASE("Sum multiple samples", "[unit][diagnostics]") {
  std::vector<std::vector<double>> samples = {
      {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {2.0, 3.0, 4.0}};

  auto result = compute_average(samples, OutputAvgType::SUM);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(6.0));
  REQUIRE(result[1] == Approx(9.0));
  REQUIRE(result[2] == Approx(12.0));
}

TEST_CASE("Std constant values", "[unit][diagnostics]") {
  // When all values are the same, STD should be 0
  std::vector<std::vector<double>> samples = {
      {5.0, 5.0, 5.0}, {5.0, 5.0, 5.0}, {5.0, 5.0, 5.0}};

  auto result = compute_average(samples, OutputAvgType::STD);

  REQUIRE(result.size() == 3);
  REQUIRE(result[0] == Approx(0.0).margin(1e-10));
  REQUIRE(result[1] == Approx(0.0).margin(1e-10));
  REQUIRE(result[2] == Approx(0.0).margin(1e-10));
}

TEST_CASE("Std variable values", "[unit][diagnostics]") {
  // For values [1, 2, 3], mean = 2, variance = (1+0+1)/3 = 2/3
  std::vector<std::vector<double>> samples = {{1.0}, {2.0}, {3.0}};

  auto result = compute_average(samples, OutputAvgType::STD);

  REQUIRE(result.size() == 1);
  double expected_std = std::sqrt(2.0 / 3.0);
  REQUIRE(result[0] == Approx(expected_std).margin(1e-10));
}
