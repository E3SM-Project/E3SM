#ifndef EMULATOR_TEST_HPP
#define EMULATOR_TEST_HPP

#include <algorithm>
#include <cstddef>
#include <emulator.hpp>
#include <process.hpp>
#include <span>
#include <stdexcept>
#include <string_view>

namespace emulator::tests {

class ForcingEmulator : public Emulator {
public:
  ForcingEmulator(CreateConfig cfg) : Emulator(cfg) {}

  static constexpr double init_val = 1.0;
  std::span<const double> temperature() { return field_view("temperature"); }
  std::span<const double> soil_moisture() {
    return field_view("soil_moisture");
  }
  std::span<const double> surface_flux() { return field_view("surface_flux"); }

  void init_data() override {
    std::ranges::fill(field_view("temperature"), init_val);
    std::ranges::fill(field_view("soil_moisture"), init_val);
    std::ranges::fill(field_view("surface_flux"), init_val);
  }

protected:
  void run_impl(const ProcessRunOpts opts) override {

    static_cast<void>(opts);

    auto temperature = field_view("temperature");
    auto soil_moisture = field_view("soil_moisture");

    if (temperature.empty()) {
      throw std::runtime_error("TestEmulator has no temperature samples");
    }

    if (soil_moisture.size() % temperature.size() != 0) {
      throw std::runtime_error("soil_moisture size is not divisible by "
                               "the number of columns");
    }

    const std::size_t ncols = temperature.size();
    const std::size_t nlevels = soil_moisture.size() / ncols;

    for (std::size_t column = 0; column < ncols; ++column) {
      temperature[column] = 280.0 + static_cast<double>(column);

      for (std::size_t level = 0; level < nlevels; ++level) {
        const std::size_t index = column * nlevels + level;

        soil_moisture[index] = 0.10 + 0.01 * static_cast<double>(column) +
                               0.001 * static_cast<double>(level);
      }
    }
  }
};

class TestEmulator : public Emulator {
public:
  TestEmulator(CreateConfig cfg) : Emulator(cfg) {}

  static constexpr double init_val = 1.0;
  std::span<double> temperature() { return field_view("temperature"); }
  std::span<double> soil_moisture() { return field_view("soil_moisture"); }
  std::span<double> surface_flux() { return field_view("surface_flux"); }
  void init_data() override {
    std::ranges::fill(field_view("temperature"), init_val);
    std::ranges::fill(field_view("soil_moisture"), init_val);
    std::ranges::fill(field_view("surface_flux"), init_val);
  }

protected:
  void run_impl(const ProcessRunOpts opts) override {
    static_cast<void>(opts);

    const auto temperature = field_view("temperature");
    const auto soil_moisture = field_view("soil_moisture");
    auto surface_flux = field_view("surface_flux");

    if (temperature.size() != surface_flux.size()) {
      throw std::runtime_error("temperature and surface_flux must have "
                               "the same column extent");
    }

    if (temperature.empty() || soil_moisture.size() % temperature.size() != 0) {
      throw std::runtime_error("Invalid soil_moisture layout");
    }

    const std::size_t ncols = temperature.size();
    const std::size_t nlevels = soil_moisture.size() / ncols;

    for (std::size_t column = 0; column < ncols; ++column) {
      double mean_soil_moisture = 0.0;

      for (std::size_t level = 0; level < nlevels; ++level) {
        const std::size_t index = column * nlevels + level;
        mean_soil_moisture += soil_moisture[index];
      }

      mean_soil_moisture /= static_cast<double>(nlevels);

      surface_flux[column] =
          2.0 * temperature[column] + 100.0 * mean_soil_moisture;
    }
  }
};

GridDesc make_test_grid() {
  return GridDesc{
      .grid_type = 0,
      .nx = 4,
      .ny = 1,
      .dimensions =
          {
              {
                  .name = "column",
                  .local_extent = 4,
                  .global_extent = 4,
              },
              {
                  .name = "soil_levels",
                  .local_extent = 3,
                  .global_extent = 3,
              },
          },
      .col_gids = {0, 1, 2, 3},
      .lat = {10.0, 20.0, 30.0, 40.0},
      .lon = {-100.0, -99.0, -98.0, -97.0},
      .area = {1.0, 1.0, 1.0, 1.0},
  };
}

} // namespace emulator::tests
#endif
