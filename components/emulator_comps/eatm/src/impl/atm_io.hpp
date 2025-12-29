/**
 * @file atm_io.hpp
 * @brief I/O functions for atmosphere emulator initial conditions.
 *
 * Provides functions to read initial conditions from NetCDF files
 * (ERA5/ACE2 format). Restart file I/O is handled by EmulatorOutputManager.
 */

#ifndef ATM_IO_HPP
#define ATM_IO_HPP

#include "../../../common/src/emulator_logger.hpp"
#include "atm_field_manager.hpp"
#include <string>
#include <vector>

namespace emulator {
namespace impl {

/**
 * @brief Read initial conditions from a NetCDF file.
 *
 * Reads atmosphere state variables from an IC file (typically ERA5 or
 * ACE2 format) and distributes the data to local partitions based on
 * the column global IDs.
 *
 * ## Supported Variables
 * Variables are read based on the `required_vars` list. Missing variables
 * are handled with defaults:
 *
 * - `global_mean_co2`: Default 415e-6
 * - `HGTsfc`: Default 0.0
 * - `DSWRFtoa`: Default 1361.0 (solar constant)
 * - Other missing vars: Default 0.0
 *
 * @param filename Path to NetCDF IC file
 * @param num_global_cols Total global columns
 * @param num_local_cols Local columns on this rank
 * @param col_gids Global IDs of local columns (1-based)
 * @param lat Latitude of each local column [radians]
 * @param fields Field manager to populate
 * @param required_vars List of variable names to read
 * @param logger Logger for status messages
 * @param is_root True if this is rank 0 (for logging)
 * @return true if successful
 */
bool read_atm_initial_conditions(const std::string &filename,
                                 int num_global_cols, int num_local_cols,
                                 const std::vector<int> &col_gids,
                                 const std::vector<double> &lat,
                                 AtmFieldManager &fields,
                                 const std::vector<std::string> &required_vars,
                                 Logger &logger, bool is_root);

} // namespace impl
} // namespace emulator

#endif // ATM_IO_HPP
