/**
 * @file samudrace_demo.cpp
 * @brief Runs the SamudrACE atmosphere and ocean through the Harness.
 *
 * Plays the part of the coupled system: it owns the external buffers, fills them
 * with the initial state and forcing that prepare_demo.py wrote, steps the
 * emulators, and writes the outputs for compare_demo.py. The atmosphere runs every
 * base step, and the ocean every n_inner_steps base steps.
 */
#include <algorithm>
#include <mpi.h>

#include <emulator.hpp>
#include <harness.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace e3sm::emulator;

namespace {

constexpr std::size_t n_lat = 180;
constexpr std::size_t n_lon = 360;
constexpr std::size_t n_values = n_lat * n_lon;
static constexpr std::string CO2_FIELD_NAME = "global_mean_co2";
constexpr float atmosphere_dt = 21600.f;

void read_field(const fs::path& path, std::vector<double>& values) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Could not open " + path.string());
  }
  file.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(double)));
  if (file.gcount() !=
      static_cast<std::streamsize>(values.size() * sizeof(double))) {
    throw std::runtime_error("Too little data in " + path.string());
  }
}

void write_field(const fs::path& path, const std::vector<double>& values) {
  fs::create_directories(path.parent_path());
  std::ofstream file(path, std::ios::binary);
  file.write(reinterpret_cast<const char*>(values.data()),
             static_cast<std::streamsize>(values.size() * sizeof(double)));
  if (!file) {
    throw std::runtime_error("Could not write " + path.string());
  }
}

GridDesc make_grid(MPI_Comm comm) {
  GridDesc grid{};
  grid.nx = n_lon;
  grid.ny = n_lat;
  grid.num_global_columns = n_values;
  grid.num_local_columns = n_values;
  grid.comm = comm;
  grid.dimensions = {{"lat", n_lat, n_lat}, {"lon", n_lon, n_lon}};
  return grid;
}

int read_count(const fs::path& path) {
  int count = 0;
  std::ifstream(path) >> count;
  if (count < 1) {
    throw std::runtime_error("Could not read a count from " + path.string());
  }
  return count;
}

using Buffers = std::map<std::string, std::vector<double>>;

/// Fills the buffers named by the .bin files in a directory, if it exists.
void load_directory(const fs::path& directory, Buffers& buffers) {
  if (!fs::exists(directory)) {
    return;
  }
  for (const auto& entry : fs::directory_iterator(directory)) {
    read_field(entry.path(), buffers.at(entry.path().stem().string()));
  }
}

void write_outputs(const fs::path& directory, Emulator& emulator,
                   const EmulatorConfig& config) {
  for (const auto& name : config.output_order) {
    const auto values = emulator.get_field(name, ModelFieldRole::Output);
    write_field(directory / (name + ".bin"),
                std::vector<double>(values.begin(), values.end()));
  }
}

void run(const fs::path& run_dir) {
  const int n_coupled_steps = read_count(run_dir / "n_coupled_steps.txt");
  const int n_inner_steps = read_count(run_dir / "n_inner_steps.txt");

  const auto atmosphere = read_emulator_config(run_dir / "atmosphere_config.yaml");
  const auto ocean = read_emulator_config(run_dir / "ocean_config.yaml");

  // The coupled system's side of the exchange: one buffer per field. A field
  // that is both an input and an output shares a buffer, which carries its
  // state from one step to the next. The plan decides which ones are used.
  Buffers buffers;
  std::set<std::string> input_names;
  std::set<std::string> output_names;
  for (const auto* config : {&atmosphere, &ocean}) {
    input_names.insert(config->input_order.begin(), config->input_order.end());
    output_names.insert(config->output_order.begin(), config->output_order.end());
  }
  ExternalBindings bindings;
  for (const auto& name : input_names) {
    buffers.try_emplace(name, n_values, 0.0);
    bindings.from_external.push_back({name, buffers.at(name)});
  }
  for (const auto& name : output_names) {
    buffers.try_emplace(name, n_values, 0.0);
    bindings.to_external.push_back({name, buffers.at(name)});
  }
  load_directory(run_dir / "initial_state", buffers);
  load_directory(run_dir / "inputs" / "fixed", buffers);

  // Routes from the ocean to the atmosphere start from the ocean's initial state,
  // because the atmosphere runs first.
  Buffers route_storage;
  RouteBufferMap route_buffers;
  for (const auto& entry :
       fs::directory_iterator(run_dir / "route_initial_state")) {
    const auto route = entry.path().stem().string();
    auto& storage = route_storage.try_emplace(route, n_values, 0.0).first->second;
    read_field(entry.path(), storage);
    route_buffers.emplace(route, storage);
  }

  // import_to_external runs every base step, and would overwrite the ocean's state
  // with its empty _out route before the ocean first runs. Using the state's own
  // buffer for both routes means there's nothing to copy.
  for (const auto& name : ocean.output_order) {
    if (std::ranges::find(ocean.input_order, name) != ocean.input_order.end()) {
      route_buffers.emplace(name + "_out", buffers.at(name));
      route_buffers.emplace(name + "_in", buffers.at(name));
    }
  }

  Harness harness;
  harness.transfer_bindings(std::move(bindings));
  harness.add_process(ProcessDescriptor{
      .run_opts = ProcessRunOpts{.dt = atmosphere_dt, .coupling_interval = 1},
      .grid = make_grid(MPI_COMM_WORLD),
      .emulator_config = atmosphere,
  });
  harness.add_process(ProcessDescriptor{
      .run_opts = ProcessRunOpts{.dt = atmosphere_dt * n_inner_steps,
                                 .coupling_interval = n_inner_steps},
      .grid = make_grid(MPI_COMM_WORLD),
      .emulator_config = ocean,
  });
  harness.initialize(run_dir / "coupling_plan.yaml", route_buffers);

  auto& atmosphere_emulator =
      harness.get_process<Emulator>(std::string_view(atmosphere.name));
  auto& ocean_emulator =
      harness.get_process<Emulator>(std::string_view(ocean.name));

  auto co2 = [](int step)->double{return 3.4545895385742185e-4 + step * 9.8904126311004355e-10;};

  const int n_steps = n_coupled_steps * n_inner_steps;
  int ocean_step = 0;
  for (int step = 0; step < n_steps; ++step) {
    const auto step_name = "step_" + std::to_string(step);
    load_directory(run_dir / "inputs" / step_name, buffers);
    if (step % n_inner_steps == 0) {
      load_directory(run_dir / "inputs" /
                         ("coupled_step_" + std::to_string(ocean_step)),
                     buffers);
    }
    std::ranges::fill(buffers.at(CO2_FIELD_NAME),co2(step));

    harness.export_from_external("cpl7");
    if (harness.run_process(atmosphere.name)) {
      write_outputs(run_dir / "outputs" / "atmosphere" / step_name,
                    atmosphere_emulator, atmosphere);
    }
    if (harness.run_process(ocean.name)) {
      write_outputs(run_dir / "outputs" / "ocean" /
                        ("step_" + std::to_string(ocean_step)),
                    ocean_emulator, ocean);
      std::cout << "Finished coupled step " << ocean_step << std::endl;
      ++ocean_step;
    }
    harness.import_to_external("cpl7");
  }
}

} // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int status = 0;
  try {
    if (argc != 2) {
      throw std::runtime_error("Usage: samudrace_demo RUN_DIR");
    }
    run(argv[1]);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    status = 1;
  }
  MPI_Finalize();
  return status;
}
