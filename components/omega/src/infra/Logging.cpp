//===-- infra/Logging.cpp - Implements Omega Logging functions --*- C++ -*-===//
//
/// \file
/// \brief implements Omega logging functions
///
/// This implements Omega logging initialization.
//
//===----------------------------------------------------------------------===//

#include "Logging.h"
#include <iostream>
#include <spdlog/sinks/basic_file_sink.h>

namespace OMEGA {

void initLogging(std::shared_ptr<spdlog::logger> Logger) {

   try {
      Logger->flush_on(spdlog::level::warn);
      spdlog::set_default_logger(Logger);

   } catch (spdlog::spdlog_ex const &Ex) {
      std::cout << "Log init failed: " << Ex.what() << std::endl;
   }
}

void initLogging(std::string const &LogFilePath) {
   auto Logger = spdlog::basic_logger_mt("*", LogFilePath);
   initLogging(Logger);
}

} // namespace OMEGA
