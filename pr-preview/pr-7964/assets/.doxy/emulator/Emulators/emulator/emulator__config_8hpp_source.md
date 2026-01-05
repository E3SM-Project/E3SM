

# File emulator\_config.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_config.hpp**](emulator__config_8hpp.md)

[Go to the documentation of this file](emulator__config_8hpp.md)


```C++


#ifndef EMULATOR_CONFIG_HPP
#define EMULATOR_CONFIG_HPP

#include "emulator_diagnostics.hpp"
#include <string>
#include <vector>

namespace emulator {

struct BuildConfig {
  std::string grid_name;                  
  std::string inference_backend = "stub"; 
};

struct RuntimeConfig {
  std::string model_path; 
  std::string ic_file;    
  bool enabled = true;    
};

struct ModelIOConfig {
  std::vector<std::string> input_variables;  
  std::vector<std::string> output_variables; 

  bool spatial_mode = true;
};

struct CouplingConfig {
  bool zero_init_exports = true; 
  bool debug = false;            
};

struct EmulatorConfig {
  BuildConfig build;            
  RuntimeConfig runtime;        
  ModelIOConfig model_io;       
  CouplingConfig coupling;      
  DiagnosticConfig diagnostics; 

  // Legacy accessors for compatibility (deprecated)
  std::string grid_file; 
};

EmulatorConfig parse_emulator_config(const std::string &yaml_file,
                                     const std::string &section_name);

EmulatorConfig
parse_emulator_config_with_defaults(const std::string &yaml_file,
                                    const std::string &section_name,
                                    bool verbose = false);

} // namespace emulator

#endif // EMULATOR_CONFIG_HPP
```


