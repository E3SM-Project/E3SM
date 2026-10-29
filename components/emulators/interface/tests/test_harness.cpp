#include "emulator_config.hpp"
#include "inference_backend.hpp"
#include <memory>
#include <string_view>
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "harness.hpp"
#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <filesystem>
#include <process.hpp>

enum class Layout {
  Strided,
  Contiguous,
};
using namespace e3sm::emulator;
static constexpr double sentinel = -999;
static constexpr std::size_t ncols = 2;
static constexpr std::size_t nlevels = 3;
static constexpr int nsteps = 3;

static std::size_t input_size() { return ncols * (1 + nlevels); }

struct FakeExt {
  FakeExt(Layout layout);

  void set_forcing(std::span<const double> temperature,
                   std::span<const double> moisture, std::size_t step);

  void run(std::string_view);

  std::vector<double> surface_flux();

  Layout layout;
  std::vector<double> x2l;
  std::vector<double> l2x;
  Harness harness;

  std::size_t temperature_index(std::size_t column) const {
    if (layout == Layout::Contiguous) {
      return column;
    }
    return column * (1 + nlevels);
  }

  std::size_t moisture_index(std::size_t column, std::size_t level) const {
    if (layout == Layout::Contiguous) {
      return ncols + column * nlevels + level;
    }
    return column * (1 + nlevels) + 1 + level;
  }

  std::size_t flux_index(std::size_t column) const {
    if (layout == Layout::Contiguous) {
      return column;
    }
    return column * 2;
  }

  std::size_t unused_index(std::size_t column) const {
    if (layout == Layout::Contiguous) {
      return ncols + column;
    }
    return column * 2 + 1;
  }
};
void FakeExt::run(std::string_view name) { harness.run_process(name); }

void FakeExt::set_forcing(std::span<const double> temperature,
                          std::span<const double> moisture, std::size_t step) {
  for (std::size_t column = 0; column < ncols; ++column) {
    x2l[temperature_index(column)] = temperature[column] + double(step);

    for (std::size_t level = 0; level < nlevels; ++level) {
      x2l[moisture_index(column, level)] =
          moisture[column * nlevels + level] - double(step) / 10.0;
    }
  }
}

std::vector<double> FakeExt::surface_flux() {
  std::vector<double> result(ncols);

  for (std::size_t column = 0; column < ncols; ++column) {
    result[column] = l2x[flux_index(column)];
  }

  return result;
}

TEST_CASE("Harness preserves results across external buffer layouts",
          "[harness][libtorch]") {

  auto land_config =
      read_emulator_config(TEST_DATA_DIR "/land_emulator_config.yaml");

  auto p_name = land_config.name;
  FakeExt direct{Layout::Contiguous};
  // FakeExt strided{Layout::Strided};

  std::vector<double> temperature{280.0, 281.0};
  std::vector<double> moisture{
      0.10, 0.11, 0.12, 0.20, 0.21, 0.22,
  };
  auto direct_flux = direct.surface_flux();

  ExternalBindings bindings{
      .from_external =
          {
              {"temperature", std::span<double>{temperature}},
              {"soil_moisture", std::span<double>{moisture}},
          },
      .to_external =
          {
              {"surface_flux", std::span<double>{direct_flux}},
          },
  };

direct.harness.transfer_bindings(bindings)

  for (int step = 0; step < nsteps; ++step) {
    CAPTURE(step);

    direct.set_forcing(temperature, moisture, step);
    // strided.set_forcing(temperature, moisture, step);

    direct.run(p_name);
    // strided.run(p_name);

    const std::vector<double> direct_flux = direct.surface_flux();
    // const std::vector<double> staged_flux = strided.surface_flux();

    REQUIRE(direct_flux.size() == ncols);
    // REQUIRE(staged_flux.size() == ncols);

    for (std::size_t column = 0; column < ncols; ++column) {
      double mean = 0.0;
      for (std::size_t level = 0; level < nlevels; ++level) {
        mean += moisture[column * nlevels + level];
      }
      mean /= static_cast<double>(nlevels);

      const double expected = 2.0 * temperature[column] + 100.0 * mean;

      REQUIRE(direct_flux[column] == Approx(expected));
      // REQUIRE(staged_flux[column] == Approx(expected));

      // REQUIRE(direct_flux[column] == staged_flux[column]);
    }

    for (std::size_t column = 0; column < ncols; ++column) {
      temperature[column] += 0.001 * direct_flux[column];
    }
  }
}
