

# File emulator\_diagnostics.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_diagnostics.hpp**](emulator__diagnostics_8hpp.md)

[Go to the documentation of this file](emulator__diagnostics_8hpp.md)


```C++


#ifndef EMULATOR_DIAGNOSTICS_HPP
#define EMULATOR_DIAGNOSTICS_HPP

#include <string>
#include <vector>

namespace emulator {

// ============================================================================
// Enums
// ============================================================================

enum class FrequencyUnit {
  NSTEPS,  
  NSECS,   
  NMINS,   
  NHOURS,  
  NDAYS,   
  NMONTHS, 
  NYEARS,  
  NONE     
};

enum class OutputAvgType {
  INSTANT, 
  AVERAGE, 
  MIN,     
  MAX,     
  STD,     
  SUM      
};

enum class OutputPrecision {
  FLOAT32, 
  FLOAT64  
};

enum class FileType {
  HISTORY,        
  RESTART,        
  HISTORY_RESTART 
};

// ============================================================================
// Configuration Structures
// ============================================================================

struct OutputStreamConfig {
  std::string stream_name = "h0";           
  std::string filename_prefix = "emulator"; 
  std::vector<std::string> fields;          

  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;
  int frequency = 1; 

  OutputAvgType avg_type = OutputAvgType::INSTANT;
  OutputPrecision precision = OutputPrecision::FLOAT32;
  int max_snapshots_per_file = 1; 
};

struct RestartConfig {
  bool enabled = true;
  std::string filename_prefix = "emulator.atm.r";
  FrequencyUnit frequency_unit = FrequencyUnit::NDAYS;
  int frequency = 1;
};

struct HistoryRestartConfig {
  bool enabled = true;
  std::string filename_prefix = "emulator.atm.rh";
};

struct DiagnosticConfig {
  std::vector<OutputStreamConfig> history_streams; 
  RestartConfig restart;                           
  HistoryRestartConfig history_restart;            
};

// ============================================================================
// Conversion Utilities
// ============================================================================

FrequencyUnit str_to_freq_unit(const std::string &s);

std::string freq_unit_to_str(FrequencyUnit u);

OutputAvgType str_to_avg_type(const std::string &s);

std::string avg_type_to_str(OutputAvgType t);

OutputPrecision str_to_precision(const std::string &s);

std::string precision_to_str(OutputPrecision p);

std::string file_type_suffix(FileType t);

} // namespace emulator

#endif // EMULATOR_DIAGNOSTICS_HPP
```


