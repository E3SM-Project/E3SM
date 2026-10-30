#include "emulator_config.hpp"
#include <algorithm>
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
static constexpr std::size_t nfields = 3;
static constexpr int nsteps = 3;

inline GridDesc make_test_grid() {
  return GridDesc{
      .grid_type = 0,
      .nx = 4,
      .ny = 1,
      .dimensions =
          {
              {
                  .name = "column",
                  .local_extent = ncols,
                  .global_extent = ncols,
              },
              {
                  .name = "soil_levels",
                  .local_extent = nlevels,
                  .global_extent = nlevels,
              },
          },
      .col_gids = {0, 1, 2, 3},
      .lat = {10.0, 20.0, 30.0, 40.0},
      .lon = {-100.0, -99.0, -98.0, -97.0},
      .area = {1.0, 1.0, 1.0, 1.0},
  };
}

struct FieldIndex {
  std::size_t offset;
  std::size_t extent;
};
using FieldDim = std::vector<std::size_t>;
struct FakeExt {
  FakeExt(Layout layout, FieldDim x2l_dims, FieldDim l2x_dims)
      : layout(layout) {
    std::size_t offset = 0;

    indices.reserve(nfields);
    for (std::size_t d : x2l_dims) {
      indices.push_back({offset, d});
      offset += d;
    }
    x2l.resize(offset);

    offset = 0;
    for (std::size_t d : l2x_dims) {
      indices.push_back({offset, d});
      offset += d;
    }
    l2x.resize(offset);
  };

  void set_forcing(std::span<const double> temperature,
                   std::span<const double> moisture);

  void run(std::string_view);

  std::span<double> surface_flux();
  std::span<double> temperature();
  std::span<double> moisture();

  Layout layout;
  std::vector<double> x2l;
  std::vector<double> l2x;
  std::vector<FieldIndex> indices;
  Harness harness;

  FieldIndex temperature_index() const { return indices[0]; }
  FieldIndex moisture_index() const { return indices[1]; }
  FieldIndex flux_index() const { return indices[2]; }
};
void FakeExt::run(std::string_view name) { harness.run_process(name); }

void FakeExt::set_forcing(std::span<const double> temperatures,
                          std::span<const double> moistures) {
  auto temperature_view = temperature();
  auto moisture_view = moisture();

  if (temperatures.size() != temperature_view.size() ||
      moistures.size() != moisture_view.size()) {
    throw std::invalid_argument("Incorrect forcing sizes");
  }
  std::ranges::copy(temperatures, temperature_view.begin());
  std::ranges::copy(moistures, moisture_view.begin());
}

std::span<double> FakeExt::temperature() {
  return std::span<double>{x2l.data() + indices[0].offset, indices[0].extent};
}
std::span<double> FakeExt::moisture() {
  return std::span<double>{x2l.data() + indices[1].offset, indices[1].extent};
}
std::span<double> FakeExt::surface_flux() {
  return std::span<double>{l2x.data() + indices[2].offset, indices[2].extent};
}

TEST_CASE("Harness preserves results across external buffer layouts",
          "[harness][libtorch]") {

  auto land_config =
      read_emulator_config(TEST_DATA_DIR "/land_emulator_config.yaml");
  auto plan = TEST_DATA_DIR "/test_harness_plan.yaml";

  auto p_name = land_config.name;
  FieldDim x2l_dims{ncols, ncols * nlevels};
  FieldDim l2x_dims{ncols};
  FakeExt direct{Layout::Contiguous, x2l_dims, l2x_dims};
  // FakeExt strided{Layout::Strided};

  REQUIRE(direct.x2l.size() == ncols * (1 + nlevels));
  REQUIRE(direct.l2x.size() == ncols);

  std::vector<double> temperature;
  std::vector<double> moisture;
  temperature.resize(ncols);
  moisture.resize(ncols * nlevels);
  for (std::size_t i = 0; i < ncols; ++i) {
    temperature[i] = 280.0 + i;
    for (std::size_t l = 0; l < nlevels; ++l) {
      moisture[i * nlevels + l] = 0.10 * (i + 1) + l;
    }
  }

  ExternalBindings bindings{
      .to_external =
          {
              {"surface_flux", direct.surface_flux()},
          },
      .from_external =
          {
              {"temperature", direct.temperature()},
              {"soil_moisture", direct.moisture()},
          },
  };

  direct.harness.transfer_bindings(bindings);
  direct.harness.add_process(ProcessDescriptor{
      .run_opts = ProcessRunOpts{.dt = 1., .coupling_interval = 1},
      .grid = make_test_grid(),
      .emulator_config = land_config,
  });
  direct.harness.initialize(
      plan, RouteBufferMap{{"surface_flux", direct.surface_flux()},
                           {"temperature", direct.temperature()},
                           {"soil_moisture", direct.moisture()}});

  Emulator& land_emu =
      direct.harness.get_process<Emulator>(std::string_view(p_name));

  std::vector<double> eflux;
  eflux.resize(ncols);

  for (int step = 0; step < nsteps; ++step) {
    CAPTURE(step);

    direct.set_forcing(temperature, moisture);
    direct.harness.export_from_external("cpl7");
    direct.run(p_name);
    direct.harness.import_to_external("cpl7");

    const auto direct_flux = direct.surface_flux();

    REQUIRE(direct_flux.size() == ncols);

    auto land_flux = land_emu.get_field("surface_flux");
    // NOTE: This should relfect what the libtorch model actually calculates!!
    //      return 2.0 * temperature + 100.0 * soil_moisture.mean(dim=1)
    for (std::size_t column = 0; column < ncols; ++column) {
      double mean = 0.0;
      for (std::size_t level = 0; level < nlevels; ++level) {
        mean += moisture[column * nlevels + level];
      }
      mean /= static_cast<double>(nlevels);

      eflux[column] = 2.0 * temperature[column] + 100.0 * mean;
      REQUIRE(eflux[column] == land_flux[column]);
      REQUIRE(direct_flux[column] == eflux[column]);
    }
    // feeback for next step
    for (std::size_t column = 0; column < ncols; ++column) {
      temperature[column] += 0.001 * direct_flux[column];
    }
  }
}
