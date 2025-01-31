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

#if defined(OMEGA_LOG_LEVEL)
#define SPDLOG_ACTIVE_LEVEL OMEGA_LOG_LEVEL
#else
// default level: info(2)
#define SPDLOG_ACTIVE_LEVEL 2
#endif

#include "DataTypes.h"
#include <spdlog/spdlog.h>
#include <string>

#include "LogFormatters.h"

#if defined(OMEGA_LOG_FLUSH)
#define _LOG_FLUSH \
   spdlog::apply_all([](std::shared_ptr<spdlog::logger> l) { l->flush(); })
#else
#define _LOG_FLUSH (void)0
#endif

#define LOG_TRACE(msg, ...)                                                   \
   spdlog::trace(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOG_DEBUG(msg, ...)                                                   \
   spdlog::debug(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOG_INFO(msg, ...)                                                   \
   spdlog::info(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOG_WARN(msg, ...)                                                   \
   spdlog::warn(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOG_ERROR(msg, ...)                                                   \
   spdlog::error(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOG_CRITICAL(msg, ...)                                   \
   spdlog::critical(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), \
                    ##__VA_ARGS__);                              \
   _LOG_FLUSH

#define LOGGER_TRACE(logger, msg, ...)                                        \
   logger->trace(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOGGER_DEBUG(logger, msg, ...)                                        \
   logger->debug(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOGGER_INFO(logger, msg, ...)                                        \
   logger->info(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOGGER_WARN(logger, msg, ...)                                        \
   logger->warn(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOGGER_ERROR(logger, msg, ...)                                        \
   logger->error(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), ##__VA_ARGS__); \
   _LOG_FLUSH
#define LOGGER_CRITICAL(logger, msg, ...)                        \
   logger->critical(OMEGA::_PackLogMsg(__FILE__, __LINE__, msg), \
                    ##__VA_ARGS__);                              \
   _LOG_FLUSH

#ifndef OMEGA_LOG_TASKS
#define OMEGA_LOG_TASKS 0
#endif

namespace OMEGA {

class MachEnv;

const std::string OmegaDefaultLogfile = "omega.log";

int initLogging(const OMEGA::MachEnv *DefEnv,
                std::shared_ptr<spdlog::logger> Logger);

int initLogging(const OMEGA::MachEnv *DefEnv,
                std::string const &LogFilePath = OmegaDefaultLogfile);

std::string _PackLogMsg(const char *file, int line, const std::string &msg);

} // namespace OMEGA

#endif // OMEGA_LOG_H
