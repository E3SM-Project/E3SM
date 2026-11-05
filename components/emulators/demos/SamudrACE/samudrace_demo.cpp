/**
 * @file samudrace_demo.cpp
 * @brief Runs the SamudrACE atmosphere through the Harness.
 *
 * Plays the part of the coupled system: it owns the external buffers, fills them
 * with the initial state and forcing that prepare_demo.py wrote, steps the
 * emulator, and writes the outputs for compare_demo.py.
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

void run(const fs::path& run_dir) {
  int n_steps = 0;
  std::ifstream(run_dir / "n_steps.txt") >> n_steps;
  if (n_steps < 1) {
    throw std::runtime_error("Could not read n_steps.txt in " + run_dir.string());
  }

  auto config = read_emulator_config(run_dir / "emulator_config.yaml");

  // The coupled system's side of the exchange: one buffer per field. A field
  // that is both an input and an output shares a buffer, which carries its
  // state from one step to the next.
  std::map<std::string, std::vector<double>> buffers;

  buffers.try_emplace(CO2_FIELD_NAME,n_values,0.0);

  ExternalBindings bindings;
  for (const auto& name : config.input_order) {
    buffers.try_emplace(name, n_values, 0.0);
    bindings.from_external.push_back({name, buffers.at(name)});
  }
  for (const auto& name : config.output_order) {
    buffers.try_emplace(name, n_values, 0.0);
    bindings.to_external.push_back({name, buffers.at(name)});
  }

  std::set<std::string> forcing;
  for (const auto& name : config.input_order) {
    const auto state = run_dir / "initial_state" / (name + ".bin");
    const auto fixed = run_dir / "inputs" / (name + ".bin");
    if (fs::exists(state)) {
      read_field(state, buffers.at(name));
    } else if (fs::exists(fixed)) {
      read_field(fixed, buffers.at(name));
    } else {
      forcing.insert(name);
      std::cout << "inserting Forcing: " << name <<std::endl;
    }
  }

  Harness harness;
  harness.transfer_bindings(std::move(bindings));
  harness.add_process(ProcessDescriptor{
      .run_opts = ProcessRunOpts{.dt = 21600.f, .coupling_interval = 1},
      .grid = make_grid(MPI_COMM_WORLD),
      .emulator_config = config,
  });
  harness.initialize(run_dir / "coupling_plan.yaml");

  auto co2 = [](int step)->double{return 3.4545895385742185e-4 + step * 9.8904126311004355e-10;};

  for (int step = 0; step < n_steps; ++step) {
    const auto step_name = "step_" + std::to_string(step);
    for (const auto& name : forcing) {
      read_field(run_dir / "inputs" / step_name / (name + ".bin"),
                 buffers.at(name));
    }
    std::ranges::fill(buffers.at(CO2_FIELD_NAME),co2(step));

    harness.export_from_external("cpl7");
    harness.run_process(config.name);
    harness.import_to_external("cpl7");

    for (const auto& name : config.output_order) {
      write_field(run_dir / "outputs" / step_name / (name + ".bin"),
                  buffers.at(name));
    }
    std::cout << "Finished " << step_name << std::endl;
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
