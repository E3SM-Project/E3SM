

# File emulator\_logger.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_logger.cpp**](emulator__logger_8cpp.md)

[Go to the documentation of this file](emulator__logger_8cpp.md)


```C++


#include "emulator_logger.hpp"
#include <ctime>
#include <iomanip>
#include <sstream>

namespace emulator {

Logger::~Logger() {
  if (m_file.is_open()) {
    m_file.close();
  }
}

void Logger::set_file(const std::string &filename) {
  if (m_file.is_open()) {
    m_file.close();
  }
  if (!filename.empty()) {
    m_file.open(filename, std::ios::out | std::ios::app);
  }
}

void Logger::log(LogLevel level, const std::string &message) {
  std::stringstream ss;
  ss << "[" << get_timestamp() << "] "
     << "[" << level_to_string(level) << "] " << message;

  if (m_file.is_open()) {
    m_file << ss.str() << std::endl;
  } else {
    std::cout << ss.str() << std::endl;
  }
}

void Logger::debug(const std::string &message) {
  log(LogLevel::DEBUG, message);
}

void Logger::verbose(const std::string &message) {
  log(LogLevel::VERBOSE, message);
}

void Logger::info(const std::string &message) { log(LogLevel::INFO, message); }

void Logger::warn(const std::string &message) {
  log(LogLevel::WARNING, message);
}

void Logger::error(const std::string &message) {
  log(LogLevel::ERROR, message);
}

std::string Logger::get_timestamp() {
  auto now = std::time(nullptr);
  auto tm = *std::localtime(&now);
  std::stringstream ss;
  ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

std::string Logger::level_to_string(LogLevel level) {
  switch (level) {
  case LogLevel::DEBUG:
    return "DEBUG";
  case LogLevel::VERBOSE:
    return "VERBOSE";
  case LogLevel::INFO:
    return "INFO";
  case LogLevel::WARNING:
    return "WARNING";
  case LogLevel::ERROR:
    return "ERROR";
  default:
    return "UNKNOWN";
  }
}

} // namespace emulator
```


