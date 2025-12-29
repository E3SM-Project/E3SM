/**
 * @file emulator_comp.cpp
 * @brief Implementation of the EmulatorComp base class.
 *
 * Provides core functionality for E3SM emulator components including
 * MPI initialization, grid management, coupling setup, and time stepping.
 *
 * @see emulator_comp.hpp for class documentation
 */

#include "emulator_comp.hpp"
#include "emulator_config.hpp"
#include "emulator_io.hpp"
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace emulator {

EmulatorComp::EmulatorComp(CompType type)
    : m_type(type), m_comm(MPI_COMM_NULL), m_comp_id(-1), m_rank(-1),
      m_nprocs(0), m_run_type(0) {}

/**
 * @brief Initialize the component instance with MPI and metadata.
 */
void EmulatorComp::create_instance(MPI_Comm comm, int comp_id,
                                   const char *input_file, int run_type,
                                   int start_ymd, int start_tod) {
  m_comm = comm;
  m_comp_id = comp_id;
  m_run_type = run_type;

  // Check if MPI communicator is valid (not null pointer or MPI_COMM_NULL)
  if (comm != MPI_COMM_NULL && comm != 0) {
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_nprocs);
  } else {
    // Fallback for non-MPI runs or invalid comm
    m_rank = 0;
    m_nprocs = 1;
  }

  if (input_file) {
    m_input_file = input_file;
  }

  if (is_root()) {
    std::string msg = "[EmulatorComp] Creating component " +
                      get_comp_name(m_type) + " based on " + m_input_file;
    m_logger.info(msg);
  }

  m_current_ymd = start_ymd;
  m_current_tod = start_tod;
  m_step_count = 0;

  // Initialize PIO subsystem
  auto comp_name = get_comp_name(m_type);
  EmulatorIO::initialize(comm, comp_name);

  // Read grid information from configuration file
  if (!m_input_file.empty()) {
    read_grid_file(m_input_file);
  }
}

/**
 * @brief Read grid information from the configuration and grid files.
 *
 * Parses the YAML configuration to find the grid file path, then reads
 * grid dimensions and coordinates using PIO.
 */
void EmulatorComp::read_grid_file(const std::string &input_file) {
  auto config =
      parse_emulator_config_with_defaults(input_file, "eatm", is_root());

  if (config.grid_file.empty()) {
    if (is_root()) {
      m_logger.info(
          "[EmulatorComp] No grid_file in config, using default grid");
    }
    setup_default_grid();
    return;
  }

  // Open grid file using PIO
  int ncid = EmulatorIO::open_file(config.grid_file);
  if (ncid < 0) {
    if (is_root()) {
      m_logger.error(
          "[EmulatorComp] Failed to open grid file, using default grid");
    }
    setup_default_grid();
    return;
  }

  // Read grid_dims variable
  // TODO: Support rank-1 grids (neXpg2) in addition to rank-2 grids
  int grid_dims[2] = {0, 0};
  if (!EmulatorIO::read_var_1d_int(ncid, "grid_dims", grid_dims, 2)) {
    if (is_root()) {
      m_logger.error("[EmulatorComp] grid_dims not found");
    }
    EmulatorIO::close_file(ncid);
    setup_default_grid();
    return;
  }

  m_nx = grid_dims[0];
  m_ny = grid_dims[1];
  m_num_global_cols = m_nx * m_ny;

  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) + " grid: " +
                  std::to_string(m_nx) + " x " + std::to_string(m_ny) + " = " +
                  std::to_string(m_num_global_cols) + " columns");
  }

  // Read coordinate variables
  std::vector<double> lon_global(m_num_global_cols);
  std::vector<double> lat_global(m_num_global_cols);
  std::vector<double> area_global(m_num_global_cols);

  EmulatorIO::read_var_1d(ncid, "grid_center_lon", lon_global.data(),
                          m_num_global_cols);
  EmulatorIO::read_var_1d(ncid, "grid_center_lat", lat_global.data(),
                          m_num_global_cols);

  if (!EmulatorIO::read_var_1d(ncid, "grid_area", area_global.data(),
                               m_num_global_cols)) {
    std::fill(area_global.begin(), area_global.end(), 1.0);
  }

  EmulatorIO::close_file(ncid);

  // Distribute to local partitions
  distribute_grid_data(lon_global, lat_global, area_global);
}

/**
 * @brief Distribute global grid data to local partitions.
 *
 * Performs a simple 1D decomposition across MPI ranks, assigning
 * contiguous blocks of columns to each rank.
 */
void EmulatorComp::distribute_grid_data(
    const std::vector<double> &lon_global,
    const std::vector<double> &lat_global,
    const std::vector<double> &area_global) {
  // Simple 1D block decomposition
  int base_cols = m_num_global_cols / m_nprocs;
  int remainder = m_num_global_cols % m_nprocs;

  if (m_rank < remainder) {
    m_num_local_cols = base_cols + 1;
  } else {
    m_num_local_cols = base_cols;
  }

  int start_idx = 0;
  for (int r = 0; r < m_rank; ++r) {
    start_idx += (r < remainder) ? (base_cols + 1) : base_cols;
  }

  // Allocate and fill local arrays
  m_col_gids.resize(m_num_local_cols);
  m_lat.resize(m_num_local_cols);
  m_lon.resize(m_num_local_cols);
  m_area.resize(m_num_local_cols);

  for (int i = 0; i < m_num_local_cols; ++i) {
    int global_idx = start_idx + i;
    m_col_gids[i] = global_idx + 1; // MCT uses 1-based indexing
    m_lon[i] = lon_global[global_idx];
    m_lat[i] = lat_global[global_idx];
    m_area[i] = area_global[global_idx];
  }

  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) + " 0th rank has " +
                  std::to_string(m_num_local_cols) + " local columns");
  }
}

/**
 * @brief Set up a default 360x180 global grid.
 *
 * Used when no grid file is available or readable.
 */
void EmulatorComp::setup_default_grid() {
  m_nx = 360;
  m_ny = 180;
  m_num_global_cols = m_nx * m_ny;

  int base_cols = m_num_global_cols / m_nprocs;
  int remainder = m_num_global_cols % m_nprocs;

  m_num_local_cols = (m_rank < remainder) ? (base_cols + 1) : base_cols;

  int start_idx = 0;
  for (int r = 0; r < m_rank; ++r) {
    start_idx += (r < remainder) ? (base_cols + 1) : base_cols;
  }

  m_col_gids.resize(m_num_local_cols);
  m_lat.resize(m_num_local_cols);
  m_lon.resize(m_num_local_cols);
  m_area.resize(m_num_local_cols);

  for (int i = 0; i < m_num_local_cols; ++i) {
    int global_idx = start_idx + i;
    m_col_gids[i] = global_idx + 1;

    int ix = global_idx % m_nx;
    int iy = global_idx / m_nx;

    m_lon[i] = (ix + 0.5) * 360.0 / m_nx;
    m_lat[i] = -90.0 + (iy + 0.5) * 180.0 / m_ny;
    m_area[i] = 1.0;
  }

  if (is_root()) {
    m_logger.info("[EmulatorComp] Using default grid: " + std::to_string(m_nx) +
                  " x " + std::to_string(m_ny));
  }
}

/**
 * @brief Override grid data with values from Fortran driver.
 *
 * Used when the driver provides grid decomposition directly instead
 * of reading from a file.
 */
void EmulatorComp::set_grid_data(int nx, int ny, int num_local_cols,
                                 int num_global_cols, const int *col_gids,
                                 const double *lat, const double *lon,
                                 const double *area) {
  m_nx = nx;
  m_ny = ny;
  m_num_local_cols = num_local_cols;
  m_num_global_cols = num_global_cols;

  m_col_gids.resize(num_local_cols);
  m_lat.resize(num_local_cols);
  m_lon.resize(num_local_cols);
  m_area.resize(num_local_cols);

  for (int i = 0; i < num_local_cols; ++i) {
    m_col_gids[i] = col_gids[i];
    m_lat[i] = lat[i];
    m_lon[i] = lon[i];
    m_area[i] = area[i];
  }

  if (is_root()) {
    m_logger.info("[EmulatorComp] Grid set via Fortran: " +
                  std::to_string(m_nx) + " x " + std::to_string(m_ny));
  }
}

/**
 * @brief Set up coupling buffer pointers from the MCT layer.
 */
void EmulatorComp::setup_coupling(double *import_data, double *export_data,
                                  int num_imports, int num_exports,
                                  int field_size) {
  m_import_data = import_data;
  m_export_data = export_data;
  m_num_imports = num_imports;
  m_num_exports = num_exports;
  m_field_size = field_size;

  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) +
                  " coupling setup: " + std::to_string(num_imports) +
                  " imports, " + std::to_string(num_exports) + " exports, " +
                  "field_size=" + std::to_string(field_size));
  }
}

/**
 * @brief Initialize the component by calling the derived class implementation.
 */
void EmulatorComp::initialize() {
  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) +
                  " Initializing component...");
  }

  init_impl();
  m_initialized = true;

  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) +
                  " initialization complete.");
  }
}

/**
 * @brief Execute one time step: import → run → export.
 */
void EmulatorComp::run(int dt) {
  if (!m_initialized) {
    throw std::runtime_error("EmulatorComp::run() called before initialize()");
  }

  advance_time(dt);
  m_step_count++;

  if (is_root()) {
    // Format timestamp as YYYY-MM-DD HH:MM:SS
    int year = m_current_ymd / 10000;
    int month = (m_current_ymd % 10000) / 100;
    int day = m_current_ymd % 100;

    int hh = m_current_tod / 3600;
    int mm = (m_current_tod % 3600) / 60;
    int ss = m_current_tod % 60;

    std::stringstream ts;
    ts << std::setfill('0') << std::setw(4) << year << "-" << std::setw(2)
       << month << "-" << std::setw(2) << day << " " << std::setw(2) << hh
       << ":" << std::setw(2) << mm << ":" << std::setw(2) << ss;

    std::string runner_name = get_comp_name(m_type);

    m_logger.info("[EmulatorComp] " + runner_name + " run step " +
                  std::to_string(m_step_count) + " (time: " + ts.str() + ")");
  }

  import_from_coupler();
  run_impl(dt);
  export_to_coupler();
}

/**
 * @brief Finalize the component and release resources.
 */
void EmulatorComp::finalize() {
  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) +
                  " finalizing component...");
  }

  final_impl();
  EmulatorIO::finalize();
  m_initialized = false;

  if (is_root()) {
    m_logger.info("[EmulatorComp] " + get_comp_name(m_type) +
                  " finalization complete.");
  }
}

void EmulatorComp::get_local_col_gids(int *gids) const {
  if (gids && !m_col_gids.empty()) {
    std::copy(m_col_gids.begin(), m_col_gids.end(), gids);
  }
}

void EmulatorComp::get_cols_latlon(double *lat, double *lon) const {
  if (lat && !m_lat.empty()) {
    std::copy(m_lat.begin(), m_lat.end(), lat);
  }
  if (lon && !m_lon.empty()) {
    std::copy(m_lon.begin(), m_lon.end(), lon);
  }
}

void EmulatorComp::get_cols_area(double *area) const {
  if (area && !m_area.empty()) {
    std::copy(m_area.begin(), m_area.end(), area);
  }
}

/**
 * @brief Advance simulation time by dt seconds.
 *
 * Handles day rollover using a NO_LEAP calendar.
 */
void EmulatorComp::advance_time(int dt) {
  m_current_tod += dt;
  const int seconds_per_day = 86400;

  while (m_current_tod >= seconds_per_day) {
    m_current_tod -= seconds_per_day;

    // Advance day using NO_LEAP calendar
    int year = m_current_ymd / 10000;
    int month = (m_current_ymd % 10000) / 100;
    int day = m_current_ymd % 100;

    day++;

    // Days per month (NO_LEAP calendar, index 0 unused)
    static const int days_in_month[] = {0,  31, 28, 31, 30, 31, 30,
                                        31, 31, 30, 31, 30, 31};

    if (day > days_in_month[month]) {
      day = 1;
      month++;
      if (month > 12) {
        month = 1;
        year++;
      }
    }

    m_current_ymd = year * 10000 + month * 100 + day;
  }
}

} // namespace emulator
