

# File emulator\_output\_manager.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_output\_manager.cpp**](emulator__output__manager_8cpp.md)

[Go to the documentation of this file](emulator__output__manager_8cpp.md)


```C++


#include "emulator_output_manager.hpp"
#include "emulator_io.hpp"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace emulator {

void EmulatorOutputManager::initialize(
    const DiagnosticConfig &config, MPI_Comm comm,
    const std::vector<int> &col_gids, int nlat, int nlon,
    const std::string &case_name, const std::string &run_dir, Logger &logger) {
  m_config = config;
  m_comm = comm;
  m_col_gids = col_gids;
  m_nlat = nlat;
  m_nlon = nlon;
  m_case_name = case_name;
  m_run_dir = run_dir;
  m_logger = &logger;

  // Check if MPI communicator is valid
  if (comm != MPI_COMM_NULL) {
    MPI_Comm_rank(comm, &m_rank);
    m_is_root = (m_rank == 0);
  } else {
    // Fallback for non-MPI runs or uninitialized comm
    m_rank = 0;
    m_is_root = true;
  }

  // Initialize history streams
  for (const auto &stream_config : config.history_streams) {
    auto stream = std::make_unique<EmulatorOutputStream>();
    stream->initialize(stream_config, comm, col_gids, nlat, nlon, logger);
    m_history_streams.push_back(std::move(stream));
  }

  // Set up restart control from YAML config
  // Note: CIME REST_N/REST_OPTION are applied in buildnml, so the config
  // we receive already has CIME-synchronized restart frequency
  if (config.restart.enabled) {
    m_restart_control.frequency = config.restart.frequency;
    m_restart_control.frequency_unit = config.restart.frequency_unit;
  } else {
    m_restart_control.frequency_unit = FrequencyUnit::NONE;
  }

  m_initialized = true;

  if (m_is_root) {
    m_logger->info("Initialized output manager with " +
                   std::to_string(m_history_streams.size()) +
                   " history stream(s)");
  }
}

void EmulatorOutputManager::setup(const FieldDataProvider &fields) {
  // Get default restart fields if not specified
  if (m_restart_fields.empty()) {
    m_restart_fields = fields.get_field_names();
  }

  if (m_is_root) {
    m_logger->info("Output manager setup complete, " +
                   std::to_string(m_restart_fields.size()) +
                   " fields available for restart");
  }
}

void EmulatorOutputManager::init_timestep(int current_step, double dt) {
  // Update each history stream
  for (auto &stream : m_history_streams) {
    stream->init_timestep(current_step, dt);
  }

  // Update restart control
  if (m_restart_control.output_enabled()) {
    m_restart_control.current_step = current_step;
    if (m_restart_control.last_write_step < 0) {
      m_restart_control.compute_next_write_step(0, dt);
      m_restart_control.last_write_step = 0;
    }
  }
}

void EmulatorOutputManager::run(int current_step,
                                const FieldDataProvider &fields) {
  if (!m_initialized) {
    return;
  }

  // Run each history stream
  for (auto &stream : m_history_streams) {
    stream->run(current_step, fields, m_case_name);
  }
}

void EmulatorOutputManager::finalize() {
  // Finalize all history streams
  for (auto &stream : m_history_streams) {
    stream->finalize();
  }

  if (m_is_root && m_initialized) {
    m_logger->info("Output manager finalized");
  }
}

bool EmulatorOutputManager::is_restart_step(int step) const {
  return m_restart_control.is_write_step(step);
}

std::string
EmulatorOutputManager::generate_restart_filename(int step,
                                                 FileType file_type) const {
  std::ostringstream oss;
  oss << m_case_name << file_type_suffix(file_type) << std::setfill('0')
      << std::setw(10) << step << ".nc";
  return oss.str();
}

bool EmulatorOutputManager::write_restart(const FieldDataProvider &fields,
                                          int step) {
  if (!m_config.restart.enabled) {
    return true;
  }

  std::string filename =
      m_run_dir + "/" + generate_restart_filename(step, FileType::RESTART);

  int ncid = EmulatorIO::create_file(filename);
  if (ncid < 0) {
    if (m_is_root) {
      m_logger->error("Failed to create restart file: " + filename);
    }
    return false;
  }

  // Define lat and lon dimensions for 2D gridded output
  int lat_dimid = EmulatorIO::define_dim(ncid, "lat", m_nlat);
  int lon_dimid = EmulatorIO::define_dim(ncid, "lon", m_nlon);

  if (lat_dimid < 0 || lon_dimid < 0) {
    if (m_is_root) {
      m_logger->error("Failed to define lat/lon dimensions in restart file");
    }
    EmulatorIO::close_file(ncid);
    return false;
  }

  // Define lat and lon coordinate variables (1D)
  // NC_DOUBLE = 6 in netcdf.h
  std::vector<int> lat_dims = {lat_dimid};
  std::vector<int> lon_dims = {lon_dimid};
  EmulatorIO::define_var(ncid, "lat", 6, lat_dims);
  EmulatorIO::define_var(ncid, "lon", 6, lon_dims);

  // Define all restart field variables as 2D (lat, lon)
  std::vector<int> field_dimids = {lat_dimid, lon_dimid};
  for (const auto &field_name : m_restart_fields) {
    const auto *data = fields.get_field(field_name);
    if (data && !data->empty()) {
      EmulatorIO::define_var(ncid, field_name, 6, field_dimids);
    }
  }

  // End define mode - CRITICAL: must be done before any writes
  if (!EmulatorIO::end_def(ncid)) {
    if (m_is_root) {
      m_logger->error("Failed to end define mode for restart file");
    }
    EmulatorIO::close_file(ncid);
    return false;
  }

  // Write lat/lon coordinate values
  std::vector<double> lat_vals(m_nlat);
  double dlat = 180.0 / m_nlat;
  for (int i = 0; i < m_nlat; ++i) {
    lat_vals[i] = -90.0 + dlat * (i + 0.5);
  }
  EmulatorIO::write_var_1d(ncid, "lat", lat_vals.data(), m_nlat);

  std::vector<double> lon_vals(m_nlon);
  double dlon = 360.0 / m_nlon;
  for (int i = 0; i < m_nlon; ++i) {
    lon_vals[i] = dlon * (i + 0.5);
  }
  EmulatorIO::write_var_1d(ncid, "lon", lon_vals.data(), m_nlon);

  // Write all restart fields as 2D (lat, lon)
  for (const auto &field_name : m_restart_fields) {
    const auto *data = fields.get_field(field_name);
    if (data && !data->empty()) {
      EmulatorIO::write_var_2d(ncid, field_name, data->data(), m_nlon, m_nlat);
    }
  }

  // Write step as attribute
  // TODO: Add attribute writing to EmulatorIO

  EmulatorIO::close_file(ncid);

  // Update rpointer
  update_rpointer(filename, FileType::RESTART);

  if (m_is_root) {
    m_logger->info("Wrote restart file: " + filename + " (" +
                   std::to_string(m_nlat) + "x" + std::to_string(m_nlon) +
                   " grid)");
  }

  // Update control
  m_restart_control.last_write_step = step;
  m_restart_control.compute_next_write_step(step, m_restart_control.dt);

  return true;
}

bool EmulatorOutputManager::read_restart(const std::string &filename) {
  // Note: This modifies fields, but we need a mutable FieldDataProvider
  // For now, return false as placeholder - callers should use AtmIO directly
  if (m_is_root) {
    m_logger->info("Would read restart from: " + filename);
  }
  return false;
}

bool EmulatorOutputManager::write_history_restart(int step) {
  if (!m_config.history_restart.enabled) {
    return true;
  }

  bool success = true;

  for (size_t i = 0; i < m_history_streams.size(); ++i) {
    auto &stream = m_history_streams[i];

    if (!stream->needs_history_restart()) {
      continue;
    }

    std::ostringstream oss;
    oss << m_run_dir << "/" << m_case_name << ".atm.rh" << i << "."
        << std::setfill('0') << std::setw(10) << step << ".nc";
    std::string filename = oss.str();

    if (!stream->write_history_restart(filename)) {
      success = false;
      if (m_is_root) {
        m_logger->error("Failed to write history restart: " + filename);
      }
    }
  }

  if (success && m_is_root) {
    m_logger->info("Wrote history restart files at step " +
                   std::to_string(step));
  }

  return success;
}

bool EmulatorOutputManager::read_history_restart(const std::string &filename) {
  // Find stream index from filename and read
  // For now, placeholder
  if (m_is_root) {
    m_logger->info("Would read history restart from: " + filename);
  }
  return false;
}

void EmulatorOutputManager::update_rpointer(const std::string &restart_file,
                                            FileType file_type) {
  if (!m_is_root) {
    return;
  }

  std::string rpointer_path = m_run_dir + "/rpointer.atm";

  // Read existing rpointer if it exists
  std::vector<std::string> lines;
  {
    std::ifstream ifs(rpointer_path);
    if (ifs.is_open()) {
      std::string line;
      while (std::getline(ifs, line)) {
        // Keep lines that aren't of the same file type
        if ((file_type == FileType::RESTART &&
             line.find(".r.") == std::string::npos) ||
            (file_type == FileType::HISTORY_RESTART &&
             line.find(".rh") == std::string::npos)) {
          lines.push_back(line);
        }
      }
    }
  }

  // Add new restart file (just the basename)
  size_t pos = restart_file.rfind('/');
  std::string basename =
      (pos != std::string::npos) ? restart_file.substr(pos + 1) : restart_file;
  lines.push_back(basename);

  // Write updated rpointer
  std::ofstream ofs(rpointer_path);
  for (const auto &line : lines) {
    ofs << line << "\n";
  }

  m_logger->info("Updated rpointer.atm with " + basename);
}

std::string
EmulatorOutputManager::find_restart_file(const std::string &rpointer_dir,
                                         FileType file_type) const {
  std::string rpointer_path = rpointer_dir + "/rpointer.atm";

  std::ifstream ifs(rpointer_path);
  if (!ifs.is_open()) {
    return "";
  }

  std::string suffix = file_type_suffix(file_type);
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.find(suffix) != std::string::npos) {
      return rpointer_dir + "/" + line;
    }
  }

  return "";
}

void EmulatorOutputManager::compute_restart_timing(int current_step,
                                                   double dt) {
  m_restart_control.compute_next_write_step(current_step, dt);
}

} // namespace emulator
```


