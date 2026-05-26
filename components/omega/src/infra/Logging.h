#ifndef OMEGA_LOG_H
#define OMEGA_LOG_H
//===-- infra/Logging.h - logging macros and custom formatters --*- C++ -*-===//
//
/// \file
/// \brief Defines logging macros and spdlog custom formatters
///
/// This header defines macros for logging. In addition, it includes
/// Omega-specific log formatters.
//
//===----------------------------------------------------------------------===//

/// \def OMEGA_LOG_LEVEL
/// Preprocessor variable that determines at which level log messages will be
/// written. All levels higher than the requested level will be written.
/// Default is the info level which will write all info, warning, error and
/// critical messages. The levels are:
///   TRACE (0)
///   DEBUG (1)
///   INFO (2)
///   WARN (3)
///   ERROR (4)
///   CRITICAL (5)
#if defined(OMEGA_LOG_LEVEL)
#define SPDLOG_ACTIVE_LEVEL OMEGA_LOG_LEVEL
#else
#define SPDLOG_ACTIVE_LEVEL 2
#endif

#include "DataTypes.h"
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <string>

#include "LogFormatters.h"

#if defined(OMEGA_LOG_FLUSH)
#define _LOG_FLUSH \
   spdlog::apply_all([](std::shared_ptr<spdlog::logger> l) { l->flush(); })
#else
#define _LOG_FLUSH (void)0
#endif

/// \def LOG_TRACE(msg, ...)
/// Macro that writes a Log message when the TRACE or higher logging level is
/// enabled. The trace label, together with the filename and line number from
/// which the /// macro is called is pre-pended to the msg.
#define LOG_TRACE(msg, ...)                                    \
   spdlog::trace(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOG_DEBUG(msg, ...)
/// Macro that writes a Log message when the DEBUG or higher logging level is
/// enabled. The debug label, together with the filename and line number from
/// which the macro is called is pre-pended to the msg.
#define LOG_DEBUG(msg, ...)                                    \
   spdlog::debug(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOG_INFO(msg, ...)
/// Macro that writes a Log message when the INFO or higher logging level is
/// enabled. This is typically the default logging level. Unlike the other
/// macros, the message is written as provided with no additional prefix.
/// This should be used for general informational messages.
#define LOG_INFO(msg, ...)                                    \
   spdlog::info(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOG_WARN(msg, ...)
/// Macro that writes a Log message when the WARN or higher logging level is
/// enabled. The warn label, together with the filename and line number from
/// which the macro is called is pre-pended to the msg.
#define LOG_WARN(msg, ...)                                    \
   spdlog::warn(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

// original: define LOG_WARN(msg, ...)                                                   \
   spdlog::warn(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH

/// \def LOG_ERROR(msg, ...)
/// Macro that writes a Log message when the ERROR or higher logging level is
/// enabled. The error label, together with the filename and line number from
/// which the macro is called is pre-pended to the msg.
#define LOG_ERROR(msg, ...)                                    \
   spdlog::error(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOG_CRITICAL(msg, ...)
/// Macro that writes a Log message when the CRITICAL or higher logging level is
/// enabled. The critical label, together with the filename and line number from
/// which the macro is called is pre-pended to the msg.
#define LOG_CRITICAL(msg, ...)                                    \
   spdlog::critical(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_TRACE(logger, msg, ...)
/// Macro similar to LOG_TRACE but for a custom input logger.
#define LOGGER_TRACE(logger, msg, ...)                                    \
   logger->trace(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_DEBUG(logger, msg, ...)
/// Macro similar to LOG_DEBUG but for a custom input logger.
#define LOGGER_DEBUG(logger, msg, ...)                                    \
   logger->debug(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_INFO(logger, msg, ...)
/// Macro similar to LOG_INFO but for a custom input logger.
#define LOGGER_INFO(logger, msg, ...)                                    \
   logger->info(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_WARN(logger, msg, ...)
/// Macro similar to LOG_WARN but for a custom input logger.
#define LOGGER_WARN(logger, msg, ...)                                    \
   logger->warn(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_ERROR(logger, msg, ...)
/// Macro similar to LOG_ERROR but for a custom input logger.
#define LOGGER_ERROR(logger, msg, ...)                                    \
   logger->error(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def LOGGER_CRITICAL(logger, msg, ...)
/// Macro similar to LOG_CRITICAL but for a custom input logger.
#define LOGGER_CRITICAL(logger, msg, ...)                                    \
   logger->critical(                                              \
       fmt::runtime(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg)), \
       ##__VA_ARGS__);                                            \
   _LOG_FLUSH

/// \def OMEGA_LOG_TASKS
/// A preprocessor variable that defines which tasks will be writing messages.
/// The default is to write a single log from task 0. If other tasks (or ALL)
/// are provided, log messages will be written to a corresponding log file with
/// the task number appended to the log file name.
#ifndef OMEGA_LOG_TASKS
#define OMEGA_LOG_TASKS 0
#endif

namespace OMEGA {

// To prevent some circular dependencies between MachEnv and Logging
class MachEnv;

// Public variables used for logging
const std::string OmegaDefaultLogfile = "omega.log"; ///< Default log filename
static std::ofstream LogFileStream;                  ///< Default log iostream

/// Initializes Omega logging, including setting up which tasks will log
/// and add redirection of stdout/stderr to log file
int initLogging(
    const OMEGA::MachEnv *DefEnv, ///< [in] MachEnv with MPI task info
    std::string const &LogFilePath = OmegaDefaultLogfile ///< [in] file name
);

/// Alternative Logging initialization using a pre-defined custom Logger
int initLogging(
    const OMEGA::MachEnv *DefEnv,          ///< [in] MachEnv with MPI task info
    std::shared_ptr<spdlog::logger> Logger ///< [in] Logger to use
);

/// Utility function to create a log message with prefix
std::string
_PackLogMsg(const char *file, ///< [in] file where log called (cpp __FILE__)
            int line,         ///< [in] src code line where log (cpp __LINE__)
            const std::string &msg ///< [in] message text
);

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_LOG_H
