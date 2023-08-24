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

#include "LogFormatters.h"
#include <spdlog/spdlog.h>

namespace OMEGA {

const std::string OmegaDefaultLogfile = "omega.log";

void initLogging(std::shared_ptr<spdlog::logger> Logger);
void initLogging(std::string const &LogFilePath);
} // namespace OMEGA

#define LOG_LEVEL_TRACE    SPDLOG_LEVEL_TRACE
#define LOG_LEVEL_DEBUG    SPDLOG_LEVEL_DEBUG
#define LOG_LEVEL_INFO     SPDLOG_LEVEL_INFO
#define LOG_LEVEL_WARN     SPDLOG_LEVEL_WARN
#define LOG_LEVEL_ERROR    SPDLOG_LEVEL_ERROR
#define LOG_LEVEL_CRITICAL SPDLOG_LEVEL_CRITICAL
#define LOG_LEVEL_OFF      SPDLOG_LEVEL_OFF

#if !defined(LOG_ACTIVE_LEVEL)
#define LOG_ACTIVE_LEVEL LOG_LEVEL_INFO
#endif

#if defined(LOG_UNBUFFERED_LOGGING)
#define _LOG_FLUSH(logger) logger->flush()
#else
#define _LOG_FLUSH(logger) (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_TRACE
#define LOGGER_TRACE(logger, ...)            \
   SPDLOG_LOGGER_TRACE(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_TRACE(...) LOGGER_TRACE(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_TRACE(logger, ...) (void)0
#define LOG_TRACE(...)            (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_DEBUG
#define LOGGER_DEBUG(logger, ...)            \
   SPDLOG_LOGGER_DEBUG(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_DEBUG(...) LOGGER_DEBUG(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_DEBUG(logger, ...) (void)0
#define LOG_DEBUG(...)            (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_INFO
#define LOGGER_INFO(logger, ...)            \
   SPDLOG_LOGGER_INFO(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_INFO(...) LOGGER_INFO(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_INFO(logger, ...) (void)0
#define LOG_INFO(...)            (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_WARN
#define LOGGER_WARN(logger, ...)            \
   SPDLOG_LOGGER_WARN(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_WARN(...) LOGGER_WARN(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_WARN(logger, ...) (void)0
#define LOG_WARN(...)            (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_ERROR
#define LOGGER_ERROR(logger, ...)            \
   SPDLOG_LOGGER_ERROR(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_ERROR(...) LOGGER_ERROR(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_ERROR(logger, ...) (void)0
#define LOG_ERROR(...)            (void)0
#endif

#if LOG_ACTIVE_LEVEL <= LOG_LEVEL_CRITICAL
#define LOGGER_CRITICAL(logger, ...)            \
   SPDLOG_LOGGER_CRITICAL(logger, __VA_ARGS__); \
   _LOG_FLUSH(logger)
#define LOG_CRITICAL(...) LOGGER_CRITICAL(spdlog::default_logger(), __VA_ARGS__)
#else
#define LOGGER_CRITICAL(logger, ...) (void)0
#define LOG_CRITICAL(...)            (void)0
#endif

#endif // OMEGA_LOG_H
