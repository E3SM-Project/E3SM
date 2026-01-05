

# File atm\_io.cpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_io.cpp**](atm__io_8cpp.md)

[Go to the documentation of this file](atm__io_8cpp.md)


```C++


#include "atm_io.hpp"
#include "../../../common/src/emulator_io.hpp"
#include <cmath>

namespace emulator {
namespace impl {

// Physical constants (currently unused but kept for reference)
// static constexpr double RDAIR = 287.04;   // Dry air gas constant [J/K/kg]
// static constexpr double STEBOL = 5.67e-8; // Stefan-Boltzmann [W/m2/K4]

bool read_atm_initial_conditions(const std::string &filename,
                                 int num_global_cols, int num_local_cols,
                                 const std::vector<int> &col_gids,
                                 const std::vector<double> &lat,
                                 AtmFieldManager &fields,
                                 const std::vector<std::string> &required_vars,
                                 Logger &logger, bool is_root) {
  (void)lat; // Currently unused

  if (is_root) {
    logger.info("Reading IC from: " + filename);
  }

  // Open the IC file
  int ncid = EmulatorIO::open_file(filename);
  if (is_root && ncid < 0) {
    logger.warn("Failed to open IC file: " + filename);
    return false;
  }

  // Check grid dimensions (typically 180 lat x 360 lon = 64800)
  int nlat = EmulatorIO::get_dim_size(ncid, "latitude");
  int nlon = EmulatorIO::get_dim_size(ncid, "longitude");

  if (is_root) {
    logger.info("IC file grid: " + std::to_string(nlat) + " x " +
                std::to_string(nlon) + " = " + std::to_string(nlat * nlon) +
                " vs expected " + std::to_string(num_global_cols));
  }

  if (nlat * nlon != num_global_cols) {
    logger.error("Grid mismatch! Cannot read IC.");
    EmulatorIO::close_file(ncid);
    return false;
  }

  // Buffer to store global data before distribution
  std::map<std::string, std::vector<double>> global_buffers;

  for (const auto &var_name : required_vars) {
    // Ensure field exists in manager (dynamic creation if needed)
    fields.register_dynamic_field(var_name);

    // Allocate global buffer for reading
    global_buffers[var_name] = std::vector<double>(num_global_cols, 0.0);

    // Try to read from file (3D slice at time 0)
    if (EmulatorIO::read_var_3d_slice(
            ncid, var_name, global_buffers[var_name].data(), nlon, nlat, 0)) {
      if (is_root) {
        logger.info("Read " + var_name);
      }
    } else {
      // Handle missing variables with defaults
      if (var_name == "global_mean_co2") {
        std::fill(global_buffers[var_name].begin(),
                  global_buffers[var_name].end(), 415e-6);
        if (is_root) {
          logger.info("Set default " + var_name + " = 415e-6");
        }
      } else if (var_name == "HGTsfc") {
        std::fill(global_buffers[var_name].begin(),
                  global_buffers[var_name].end(), 0.0);
        if (is_root) {
          logger.info("Set default " + var_name + " = 0.0");
        }
      } else if (var_name == "DSWRFtoa") {
        std::fill(global_buffers[var_name].begin(),
                  global_buffers[var_name].end(), 1361.0);
        if (is_root) {
          logger.info("Set default " + var_name + " = 1361.0");
        }
      } else {
        if (is_root) {
          logger.warn("Failed to read " + var_name + ", using zeros.");
        }
      }
    }
  }

  EmulatorIO::close_file(ncid);

  // Distribute to local partition based on col_gids
  for (const auto &var_name : required_vars) {
    std::vector<double> *local_ptr = fields.get_field_ptr(var_name);
    const auto &global_buf = global_buffers[var_name];

    if (local_ptr) {
      for (int i = 0; i < num_local_cols; ++i) {
        int global_idx = col_gids[i] - 1; // 1-based to 0-based
        if (global_idx >= 0 && global_idx < num_global_cols) {
          (*local_ptr)[i] = global_buf[global_idx];
        }
      }
    }
  }

  if (is_root) {
    logger.info("IC generic read complete.");
  }

  return true;
}

} // namespace impl
} // namespace emulator
```


