/**
 * @file emulator_output_stream.hpp
 * @brief Single output stream handler for diagnostic output.
 *
 * Manages a single output stream including file creation, averaging,
 * and snapshot writing at configured intervals.
 */

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

/**
 * @brief Controls output timing and tracks averaging state.
 */
struct OutputControl {
  int frequency = 1;
  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;

  int nsamples_since_last_write = 0;
  int last_write_step = -1;
  int next_write_step = 0;

  // Timestep info for computing next write
  double dt = 0.0;
  int current_step = 0;

  /**
   * @brief Check if output is enabled.
   */
  bool output_enabled() const { return frequency_unit != FrequencyUnit::NONE; }

  /**
   * @brief Check if current step is a write step.
   */
  bool is_write_step(int step) const;

  /**
   * @brief Compute next write step from frequency and current step.
   */
  void compute_next_write_step(int current_step, double dt);

  /**
   * @brief Get number of seconds per frequency unit.
   */
  double seconds_per_unit() const;
};

/**
 * @brief Interface for providing field data to output streams.
 *
 * Components implement this to provide access to their field data.
 */
class FieldDataProvider {
public:
  virtual ~FieldDataProvider() = default;

  /**
   * @brief Get pointer to field data by name.
   * @param name Field name
   * @return Pointer to data vector, or nullptr if not found
   */
  virtual const std::vector<double> *
  get_field(const std::string &name) const = 0;

  /**
   * @brief Get list of all available field names.
   */
  virtual std::vector<std::string> get_field_names() const = 0;

  /**
   * @brief Get number of local columns.
   */
  virtual int get_ncols() const = 0;

  /**
   * @brief Check if a field has vertical levels (is 3D).
   * @return Number of levels, or 1 for 2D fields
   */
  virtual int get_field_nlevs(const std::string &name) const { return 1; }
};

/**
 * @brief Manages a single diagnostic output stream.
 *
 * Handles:
 *
 * - NetCDF file creation and management
 * - Accumulation buffer for averaging (AVERAGE, MIN, MAX, STD, SUM)
 * - Write timing based on configured frequency
 * - Field stacking for sliced variables
 */
class EmulatorOutputStream {
public:
  EmulatorOutputStream() = default;
  ~EmulatorOutputStream() = default;

  /**
   * @brief Initialize the output stream.
   * @param config Stream configuration
   * @param comm MPI communicator
   * @param col_gids Global IDs of local columns
   * @param nlat Number of latitudes (for global grid)
   * @param nlon Number of longitudes
   * @param logger Logger reference
   */
  void initialize(const OutputStreamConfig &config, MPI_Comm comm,
                  const std::vector<int> &col_gids, int nlat, int nlon,
                  Logger &logger);

  /**
   * @brief Called at the start of each timestep.
   * @param current_step Current model step number
   * @param dt Timestep in seconds
   */
  void init_timestep(int current_step, double dt);

  /**
   * @brief Process fields at current timestep.
   *
   * Updates averaging buffers if not an instant stream.
   * Writes output if this is a write step.
   *
   * @param current_step Current step number
   * @param fields Field data provider
   * @param case_name Case name for filename generation
   */
  void run(int current_step, const FieldDataProvider &fields,
           const std::string &case_name);

  /**
   * @brief Finalize and close any open files.
   */
  void finalize();

  /**
   * @brief Check if this step is a write step.
   */
  bool is_write_step() const {
    return m_control.is_write_step(m_control.current_step);
  }

  /**
   * @brief Get the stream configuration.
   */
  const OutputStreamConfig &config() const { return m_config; }

  // =========================================================================
  // History restart support
  // =========================================================================

  /**
   * @brief Write averaging state to history restart file.
   * @param filename Output filename
   * @return true if successful
   */
  bool write_history_restart(const std::string &filename);

  /**
   * @brief Read averaging state from history restart file.
   * @param filename Input filename
   * @return true if successful
   */
  bool read_history_restart(const std::string &filename);

  /**
   * @brief Check if stream needs history restart (non-instant averaging).
   */
  bool needs_history_restart() const {
    return m_config.avg_type != OutputAvgType::INSTANT;
  }

private:
  // =========================================================================
  // Internal methods
  // =========================================================================

  /**
   * @brief Update averaging buffers with current field values.
   */
  void update_averaging(const FieldDataProvider &fields);

  /**
   * @brief Write output to file.
   */
  void write_output(const FieldDataProvider &fields,
                    const std::string &case_name);

  /**
   * @brief Create new output file.
   */
  void setup_file(const std::string &filename);

  /**
   * @brief Generate filename for current snapshot.
   */
  std::string generate_filename(const std::string &case_name, int step) const;

  /**
   * @brief Reset averaging buffers after write.
   */
  void reset_averaging_buffers();

  /**
   * @brief Get accumulated/averaged data for a field.
   */
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
