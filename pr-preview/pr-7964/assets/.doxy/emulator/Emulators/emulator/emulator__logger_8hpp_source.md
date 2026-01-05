

# File emulator\_logger.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_logger.hpp**](emulator__logger_8hpp.md)

[Go to the documentation of this file](emulator__logger_8hpp.md)


```C++


#ifndef EMULATOR_LOGGER_HPP
#define EMULATOR_LOGGER_HPP

#include <fstream>
#include <iostream>
#include <string>

namespace emulator {

enum class LogLevel {
  DEBUG,   
  VERBOSE, 
  INFO,    
  WARNING, 
  ERROR    
};

class Logger {
public:
  Logger() = default;
  ~Logger();

  void set_file(const std::string &filename);

  void log(LogLevel level, const std::string &message);

  void debug(const std::string &message);

  void verbose(const std::string &message);

  void info(const std::string &message);

  void warn(const std::string &message);

  void error(const std::string &message);

private:
  std::string get_timestamp();
  std::string level_to_string(LogLevel level);

  std::ofstream m_file; 
};

} // namespace emulator

#endif // EMULATOR_LOGGER_HPP
```


