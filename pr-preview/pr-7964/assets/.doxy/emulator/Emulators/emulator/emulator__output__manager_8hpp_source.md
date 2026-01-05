

# File emulator\_output\_manager.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_output\_manager.hpp**](emulator__output__manager_8hpp.md)

[Go to the documentation of this file](emulator__output__manager_8hpp.md)


```C++


#ifndef EMULATOR_OUTPUT_MANAGER_HPP
#define EMULATOR_OUTPUT_MANAGER_HPP

#include "emulator_diagnostics.hpp"
#include "emulator_logger.hpp"
#include "emulator_output_stream.hpp"
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace emulator {

class EmulatorOutputManager {
public:
  EmulatorOutputManager() = default;
  ~EmulatorOutputManager() = default;

  void initialize(const DiagnosticConfig &config, MPI_Comm comm,
                  const std::vector<int> &col_gids, int nlat, int nlon,
                  const std::string &case_name, const std::string &run_dir,
                  Logger &logger);

  void setup(const FieldDataProvider &fields);

  void init_timestep(int current_step, double dt);

  void run(int current_step, const FieldDataProvider &fields);

  void finalize();

  // =========================================================================
  // Restart methods
  // =========================================================================

  bool write_restart(const FieldDataProvider &fields, int step);

  bool read_restart(const std::string &filename);

  bool write_history_restart(int step);

  bool read_history_restart(const std::string &filename);

  bool is_restart_step(int step) const;

  std::string find_restart_file(const std::string &rpointer_dir,
                                FileType file_type) const;

  // =========================================================================
  // Accessors
  // =========================================================================

  size_t num_history_streams() const { return m_history_streams.size(); }

  bool is_initialized() const { return m_initialized; }

private:
  // =========================================================================
  // Internal methods
  // =========================================================================

  void update_rpointer(const std::string &restart_file, FileType file_type);

  std::string generate_restart_filename(int step, FileType file_type) const;

  void compute_restart_timing(int current_step, double dt);

  // =========================================================================
  // Data members
  // =========================================================================

  DiagnosticConfig m_config;

  // History output streams
  std::vector<std::unique_ptr<EmulatorOutputStream>> m_history_streams;

  // Restart control
  OutputControl m_restart_control;
  std::vector<std::string> m_restart_fields; // Fields for restart

  MPI_Comm m_comm;
  int m_rank = 0;
  bool m_is_root = false;

  std::vector<int> m_col_gids;
  int m_nlat = 0;
  int m_nlon = 0;

  std::string m_case_name;
  std::string m_run_dir;

  Logger *m_logger = nullptr;
  bool m_initialized = false;
};

} // namespace emulator

#endif // EMULATOR_OUTPUT_MANAGER_HPP
```


