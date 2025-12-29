/**
 * @file emulator_logger.hpp
 * @brief Simple logging utility for emulator components.
 *
 * Provides timestamped logging to stdout or an optional file.
 * Used by EmulatorComp and derived classes for runtime diagnostics.
 */

#ifndef EMULATOR_LOGGER_HPP
#define EMULATOR_LOGGER_HPP

#include <fstream>
#include <iostream>
#include <string>

namespace emulator {

/**
 * @brief Log level enumeration.
 */
enum class LogLevel {
  DEBUG,   ///< Detailed debugging information
  VERBOSE, ///< Verbose output for development
  INFO,    ///< General informational messages
  WARNING, ///< Warning conditions
  ERROR    ///< Error conditions
};

/**
 * @brief Simple logger with optional file output.
 *
 * By default, logs to stdout. Call set_file() to redirect output
 * to a log file. Each message is prefixed with a timestamp and
 * log level.
 *
 * ## Usage
 * ```cpp
 * Logger logger;
 * logger.set_file("component.log");
 * logger.info("Component initialized");
 * logger.warn("Configuration missing, using defaults");
 * ```
 */
class Logger {
public:
  Logger() = default;
  ~Logger();

  /**
   * @brief Set the log file path.
   *
   * If a file is already open, it is closed before opening the new file.
   * If filename is empty, logging reverts to stdout.
   *
   * @param filename Path to log file (empty string for stdout)
   */
  void set_file(const std::string &filename);

  /**
   * @brief Log a message at the specified level.
   * @param level Log level
   * @param message Message to log
   */
  void log(LogLevel level, const std::string &message);

  /** @brief Log a debug message. */
  void debug(const std::string &message);

  /** @brief Log a verbose message. */
  void verbose(const std::string &message);

  /** @brief Log an informational message. */
  void info(const std::string &message);

  /** @brief Log a warning message. */
  void warn(const std::string &message);

  /** @brief Log an error message. */
  void error(const std::string &message);

private:
  std::string get_timestamp();
  std::string level_to_string(LogLevel level);

  std::ofstream m_file; ///< Optional output file
};

} // namespace emulator

#endif // EMULATOR_LOGGER_HPP
