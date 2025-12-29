/**
 * @file emulator_diagnostics.hpp
 * @brief Configuration structures for diagnostic output.
 *
 * Defines enums and structs for configuring output streams including
 * history files (.atm.h.), restart files (.atm.r.), and history restart files
 * (.atm.rh.).
 */

#ifndef EMULATOR_DIAGNOSTICS_HPP
#define EMULATOR_DIAGNOSTICS_HPP

#include <string>
#include <vector>

namespace emulator {

// ============================================================================
// Enums
// ============================================================================

/**
 * @brief Output frequency units.
 *
 * Specifies the time unit for output frequency. Follows EAMxx conventions.
 */
enum class FrequencyUnit {
  NSTEPS,  ///< Every N model steps
  NSECS,   ///< Every N seconds
  NMINS,   ///< Every N minutes
  NHOURS,  ///< Every N hours
  NDAYS,   ///< Every N days
  NMONTHS, ///< Every N months
  NYEARS,  ///< Every N years
  NONE     ///< Disabled
};

/**
 * @brief Output averaging type.
 *
 * Specifies how to combine multiple timesteps in the output window.
 */
enum class OutputAvgType {
  INSTANT, ///< No averaging, instantaneous value
  AVERAGE, ///< Mean over output window
  MIN,     ///< Minimum over output window
  MAX,     ///< Maximum over output window
  STD,     ///< Standard deviation over output window
  SUM      ///< Accumulated sum (e.g., precipitation)
};

/**
 * @brief Output data precision.
 */
enum class OutputPrecision {
  FLOAT32, ///< Single precision (default, saves space)
  FLOAT64  ///< Double precision
};

/**
 * @brief File type indicator for restart discovery.
 */
enum class FileType {
  HISTORY,        ///< Regular output (.h. files)
  RESTART,        ///< Model state for restarting (.r. files)
  HISTORY_RESTART ///< Diagnostic state for restart (.rh. files)
};

// ============================================================================
// Configuration Structures
// ============================================================================

/**
 * @brief Configuration for a single output stream.
 *
 * Each stream writes to its own set of NetCDF files with configurable
 * output frequency, averaging, and field selection.
 */
struct OutputStreamConfig {
  std::string stream_name = "h0";           ///< Stream identifier
  std::string filename_prefix = "emulator"; ///< Output filename prefix
  std::vector<std::string> fields;          ///< Fields to output

  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;
  int frequency = 1; ///< Output every N units

  OutputAvgType avg_type = OutputAvgType::INSTANT;
  OutputPrecision precision = OutputPrecision::FLOAT32;
  int max_snapshots_per_file = 1; ///< Snapshots before new file
};

/**
 * @brief Restart output configuration.
 */
struct RestartConfig {
  bool enabled = true;
  std::string filename_prefix = "emulator.atm.r";
  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;
  int frequency = 1;
};

/**
 * @brief History restart configuration for averaging buffers.
 */
struct HistoryRestartConfig {
  bool enabled = true;
  std::string filename_prefix = "emulator.atm.rh";
};

/**
 * @brief Complete diagnostic output configuration.
 *
 * Aggregates all output stream configurations. Parsed from the
 * `diagnostics` section of the YAML configuration file.
 *
 * ## Example YAML
 * ```yaml
 * diagnostics:
 *   history:
 *     - stream_name: h0
 *       filename_prefix: "eatm.h0"
 *       frequency: 1
 *       frequency_unit: ndays
 *       averaging: instant
 *       fields: [surface_temperature, PRESsfc]
 *   restart:
 *     enabled: true
 *     frequency: 1
 *     frequency_unit: ndays
 * ```
 */
struct DiagnosticConfig {
  std::vector<OutputStreamConfig> history_streams; ///< .h. file streams
  RestartConfig restart;                           ///< .r. file config
  HistoryRestartConfig history_restart;            ///< .rh. file config
};

// ============================================================================
// Conversion Utilities
// ============================================================================

/**
 * @brief Convert string to FrequencyUnit.
 * @param s String like "nsteps", "ndays", "nmonths", etc.
 * @return Corresponding FrequencyUnit (NONE if invalid)
 */
FrequencyUnit str_to_freq_unit(const std::string &s);

/**
 * @brief Convert FrequencyUnit to string.
 */
std::string freq_unit_to_str(FrequencyUnit u);

/**
 * @brief Convert string to OutputAvgType.
 * @param s String like "instant", "average", "min", "max", "std", "sum"
 * @return Corresponding OutputAvgType
 */
OutputAvgType str_to_avg_type(const std::string &s);

/**
 * @brief Convert OutputAvgType to string.
 */
std::string avg_type_to_str(OutputAvgType t);

/**
 * @brief Convert string to OutputPrecision.
 */
OutputPrecision str_to_precision(const std::string &s);

/**
 * @brief Convert OutputPrecision to string.
 */
std::string precision_to_str(OutputPrecision p);

/**
 * @brief Get file type suffix string.
 * @return ".atm.h.", ".atm.r.", or ".atm.rh."
 */
std::string file_type_suffix(FileType t);

} // namespace emulator

#endif // EMULATOR_DIAGNOSTICS_HPP
