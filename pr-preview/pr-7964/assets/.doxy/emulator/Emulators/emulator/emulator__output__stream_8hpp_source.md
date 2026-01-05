

# File emulator\_output\_stream.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_output\_stream.hpp**](emulator__output__stream_8hpp.md)

[Go to the documentation of this file](emulator__output__stream_8hpp.md)


```C++


#ifndef EMULATOR_OUTPUT_STREAM_HPP
#define EMULATOR_OUTPUT_STREAM_HPP

#include "emulator_diagnostics.hpp"
#include "emulator_logger.hpp"
#include <functional>
#include <map>
#include <mpi.h>
#include <string>
#include <vector>

namespace emulator {

// Forward declaration for field data access
class FieldDataProvider;

struct OutputControl {
  int frequency = 1;
  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;

  int nsamples_since_last_write = 0;
  int last_write_step = -1;
  int next_write_step = 0;

  // Timestep info for computing next write
  double dt = 0.0;
  int current_step = 0;

  bool output_enabled() const { return frequency_unit != FrequencyUnit::NONE; }

  bool is_write_step(int step) const;

  void compute_next_write_step(int current_step, double dt);

  double seconds_per_unit() const;
};

class FieldDataProvider {
public:
  virtual ~FieldDataProvider() = default;

  virtual const std::vector<double> *
  get_field(const std::string &name) const = 0;

  virtual std::vector<std::string> get_field_names() const = 0;

  virtual int get_ncols() const = 0;

  virtual int get_field_nlevs(const std::string &name) const { return 1; }
};

class EmulatorOutputStream {
public:
  EmulatorOutputStream() = default;
  ~EmulatorOutputStream() = default;

  void initialize(const OutputStreamConfig &config, MPI_Comm comm,
                  const std::vector<int> &col_gids, int nlat, int nlon,
                  Logger &logger);

  void init_timestep(int current_step, double dt);

  void run(int current_step, const FieldDataProvider &fields,
           const std::string &case_name);

  void finalize();

  bool is_write_step() const {
    return m_control.is_write_step(m_control.current_step);
  }

  const OutputStreamConfig &config() const { return m_config; }

  // =========================================================================
  // History restart support
  // =========================================================================

  bool write_history_restart(const std::string &filename);

  bool read_history_restart(const std::string &filename);

  bool needs_history_restart() const {
    return m_config.avg_type != OutputAvgType::INSTANT;
  }

private:
  // =========================================================================
  // Internal methods
  // =========================================================================

  void update_averaging(const FieldDataProvider &fields);

  void write_output(const FieldDataProvider &fields,
                    const std::string &case_name);

  void setup_file(const std::string &filename);

  std::string generate_filename(const std::string &case_name, int step) const;

  void reset_averaging_buffers();

  std::vector<double> get_output_data(const std::string &field_name,
                                      const FieldDataProvider &fields) const;

  // =========================================================================
  // Data members
  // =========================================================================

  OutputStreamConfig m_config;
  OutputControl m_control;

  MPI_Comm m_comm;
  int m_rank = 0;
  bool m_is_root = false;

  std::vector<int> m_col_gids;
  int m_nlat = 0;
  int m_nlon = 0;
  int m_ncols_local = 0;

  // Averaging buffers: field_name -> accumulated data
  std::map<std::string, std::vector<double>> m_avg_buffer;
  // For STD: need sum of squares
  std::map<std::string, std::vector<double>> m_avg_buffer_sq;

  // Current file state
  int m_current_file_ncid = -1;
  int m_snapshots_in_file = 0;
  std::string m_current_filename;

  Logger *m_logger = nullptr;
  bool m_initialized = false;
};

} // namespace emulator

#endif // EMULATOR_OUTPUT_STREAM_HPP
```


