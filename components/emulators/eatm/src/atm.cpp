/**
 * @file atm.cpp
 * @brief Atmosphere emulator component implementation.
 *
 * Stub implementation — fill in details for AI/ML inference,
 * coupling, and I/O.
 */

#include "atm.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <mpi.h>

namespace emulator {

EmulatorAtm::EmulatorAtm()
    : Emulator(EmulatorType::ATM_COMP, -1, "eatm") {}

void EmulatorAtm::create_instance(int comm, int comp_id,
                                  const std::string &input_file,
                                  int run_type, int start_ymd,
                                  int start_tod) {
  m_comm = comm;
  m_id = comp_id;  // set base class ID
  m_input_file = input_file;
  m_run_type = run_type;
  (void)start_ymd;
  (void)start_tod;

  // Simple configuration parsing
  if (!input_file.empty()) {
    std::ifstream ifs(input_file);
    std::string line;
    while (std::getline(ifs, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      size_t pos = line.find(':');
      if (pos != std::string::npos) {
        std::string key = line.substr(0, pos);
        std::string val = line.substr(pos + 1);
        // trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        val.erase(0, val.find_first_not_of(" \t"));
        val.erase(val.find_last_not_of(" \t") + 1);

        if (key == "nx") {
          m_nx = std::stoi(val);
        }
        if (key == "ny") {
          m_ny = std::stoi(val);
        }
        if (key == "grid") {
          // grid name identified
        }
      }
    }
  }

  // Compute global column count from grid dimensions
  if (m_nx > 0) {
    m_num_global_cols = m_nx * std::max(1, m_ny);
  }

  // If we have a global size but no local size, create a default decomposition
  if (m_num_global_cols > 0 && m_num_local_cols == 0) {
    int rank, size;
    MPI_Comm_rank(MPI_Comm_f2c(m_comm), &rank);
    MPI_Comm_size(MPI_Comm_f2c(m_comm), &size);

    int n_per_rank = m_num_global_cols / size;
    int remainder = m_num_global_cols % size;

    int start_idx = rank * n_per_rank + std::min(rank, remainder);
    m_num_local_cols = n_per_rank + (rank < remainder ? 1 : 0);

    m_col_gids.resize(m_num_local_cols);
    for (int i = 0; i < m_num_local_cols; ++i) {
      m_col_gids[i] = start_idx + i + 1; // 1-based GIDs for MCT
    }

    m_lat.assign(m_num_local_cols, 0.0);
    m_lon.assign(m_num_local_cols, 0.0);
    m_area.assign(m_num_local_cols, 1.0);
  }
}

void EmulatorAtm::set_grid_data(int nx, int ny,
                                int num_local_cols,
                                int num_global_cols,
                                const int *col_gids,
                                const double *lat,
                                const double *lon,
                                const double *area) {
  m_nx = nx;
  m_ny = ny;
  m_num_local_cols = num_local_cols;
  m_num_global_cols = num_global_cols;

  m_col_gids.assign(col_gids, col_gids + num_local_cols);
  m_lat.assign(lat, lat + num_local_cols);
  m_lon.assign(lon, lon + num_local_cols);
  m_area.assign(area, area + num_local_cols);
}

void EmulatorAtm::init_coupling_indices(
    const std::string &export_fields,
    const std::string &import_fields) {
  // TODO: Parse colon-separated MCT field lists and populate
  // m_coupling_idx with index positions.
  (void)export_fields;
  (void)import_fields;
}

void EmulatorAtm::setup_coupling(double *import_data,
                                 double *export_data,
                                 int num_imports,
                                 int num_exports,
                                 int field_size) {
  m_import_data = import_data;
  m_export_data = export_data;
  m_num_imports = num_imports;
  m_num_exports = num_exports;
  (void)field_size;
}

void EmulatorAtm::get_local_col_gids(int *gids) const {
  std::memcpy(gids, m_col_gids.data(),
              m_col_gids.size() * sizeof(int));
}

void EmulatorAtm::get_cols_latlon(double *lat, double *lon) const {
  std::memcpy(lat, m_lat.data(),
              m_lat.size() * sizeof(double));
  std::memcpy(lon, m_lon.data(),
              m_lon.size() * sizeof(double));
}

void EmulatorAtm::get_cols_area(double *area) const {
  std::memcpy(area, m_area.data(),
              m_area.size() * sizeof(double));
}

// =========================================================================
// Lifecycle implementations
// =========================================================================

void EmulatorAtm::init_impl() {
  // TODO: Load YAML configuration from m_input_file
  // TODO: Create inference backend
  // TODO: Read initial conditions
  // TODO: Set up diagnostic output manager

  // TODO: Allocate field storage
  // TODO: Export initial values to coupler
}

void EmulatorAtm::run_impl(int dt) {
  (void)dt;

  // 1. Import fields from coupler
  import_coupling_fields();

  // 2. Prepare AI model inputs
  prepare_inputs();

  // 3. TODO: Run AI inference
  // run_inference(m_fields.net_inputs, m_fields.net_outputs);

  // 4. Process AI outputs
  process_outputs();

  // 5. TODO: Diagnostic output

  // 6. Export fields to coupler
  export_coupling_fields();
}

void EmulatorAtm::final_impl() {
  // TODO: Write final restart files
  // TODO: Finalize output manager
  // TODO: Finalize inference backend

  // TODO: Deallocate field storage
  std::cout << "eatm c++ side ... bye!" << std::endl;
}

// =========================================================================
// Coupling helpers
// =========================================================================

void EmulatorAtm::import_coupling_fields() {
  // TODO: Transfer coupler import data → internal fields
}

void EmulatorAtm::export_coupling_fields() {
  // TODO: Transfer internal fields → coupler export data
}

void EmulatorAtm::prepare_inputs() {
  // TODO: Pack field data into m_fields.net_inputs tensor
  // for inference. Handle spatial_mode vs pointwise layout.
}

void EmulatorAtm::process_outputs() {
  // TODO: Unpack m_fields.net_outputs tensor into field
  // vectors. Handle spatial_mode vs pointwise layout.
}

} // namespace emulator
