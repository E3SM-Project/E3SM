/**
 * @file emulator_output_manager.hpp
 * @brief Manages all diagnostic output streams.
 *
 * Coordinates multiple output streams, restart output, and history restart
 * files. Provides unified interface for EmulatorComp to manage I/O.
 */

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

/**
 * @brief Manages all diagnostic output for an emulator component.
 *
 * ## Responsibilities
 *
 * - Multiple history output streams (.h. files)
 * - Model restart output (.r. files)
 * - History restart output (.rh. files)
 * - rpointer.atm management for restart discovery
 *
 * ## Lifecycle
 *
 * 1. `initialize()` - Parse config, set up internals
 * 2. `setup()` - Register fields from field provider
 * 3. `init_timestep()` - Called each step before run
 * 4. `run()` - Update averaging, write if needed
 * 5. `finalize()` - Close all files
 */
class EmulatorOutputManager {
public:
  EmulatorOutputManager() = default;
  ~EmulatorOutputManager() = default;

  /**
   * @brief Initialize the output manager.
   * @param config Diagnostic configuration
   * @param comm MPI communicator
   * @param col_gids Global IDs of local columns
   * @param nlat Number of latitudes
   * @param nlon Number of longitudes
   * @param case_name Case name for filenames
   * @param run_dir Run directory for output files
   * @param logger Logger reference
   */
  void initialize(const DiagnosticConfig &config, MPI_Comm comm,
                  const std::vector<int> &col_gids, int nlat, int nlon,
                  const std::string &case_name, const std::string &run_dir,
                  Logger &logger);

  /**
   * @brief Set up output streams with field information.
   * @param fields Field data provider
   */
  void setup(const FieldDataProvider &fields);

  /**
   * @brief Called at the start of each timestep.
   * @param current_step Current model step
   * @param dt Timestep in seconds
   */
  void init_timestep(int current_step, double dt);

  /**
   * @brief Process output at current timestep.
   *
   * Updates averaging buffers and writes output if needed.
   *
   * @param current_step Current model step
   * @param fields Field data provider
   */
  void run(int current_step, const FieldDataProvider &fields);

  /**
   * @brief Finalize and close all files.
   */
  void finalize();

  // =========================================================================
  // Restart methods
  // =========================================================================

  /**
   * @brief Write model restart file.
   *
   * Writes all prognostic fields to a restart file and updates rpointer.atm.
   *
   * @param fields Field data provider
   * @param step Current step (for filename)
   * @return true if successful
   */
  bool write_restart(const FieldDataProvider &fields, int step);

  /**
   * @brief Read model restart file.
   * @param filename Restart file path
   * @param fields Field data receiver (mutable)
   * @return true if successful
   */
  bool read_restart(const std::string &filename);

  /**
   * @brief Write history restart for all non-instant streams.
   * @param step Current step (for filename)
   * @return true if successful
   */
  bool write_history_restart(int step);

  /**
   * @brief Read history restart from file.
   * @param filename History restart file path
   * @return true if successful
   */
  bool read_history_restart(const std::string &filename);

  /**
   * @brief Check if restart should be written at this step.
   */
  bool is_restart_step(int step) const;

  /**
   * @brief Find restart file from rpointer.
   * @param rpointer_dir Directory containing rpointer.atm
   * @param file_type RESTART or HISTORY_RESTART
   * @return Filename if found, empty string otherwise
   */
  std::string find_restart_file(const std::string &rpointer_dir,
                                FileType file_type) const;

  // =========================================================================
  // Accessors
  // =========================================================================

  /**
   * @brief Get number of history streams.
   */
  size_t num_history_streams() const { return m_history_streams.size(); }

  /**
   * @brief Check if output manager is initialized.
   */
  bool is_initialized() const { return m_initialized; }

private:
  // =========================================================================
  // Internal methods
  // =========================================================================

  /**
   * @brief Update rpointer.atm with new restart file.
   */
  void update_rpointer(const std::string &restart_file, FileType file_type);

  /**
   * @brief Generate restart filename.
   */
  std::string generate_restart_filename(int step, FileType file_type) const;

  /**
   * @brief Compute restart control timing.
   */
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
