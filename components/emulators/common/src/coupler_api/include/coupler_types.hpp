#ifndef COUPLER_TYPES_HPP
#define COUPLER_TYPES_HPP

#include <cstddef>
extern "C" {
/**
 * @brief Configuration parameters for creating emulator instance.
 *
 * Fields:
 *  - f_comm: MPI communicator from Fortran
 *  - comp_id: component id
 *  - run_type: cold start, restart, etc...
 *  - start_ymd: simulation start date
 *  - start_tod: time of day in seconds
 *  - start_tod: simulation start time of day in seconds
 *  - case_start_ymd: case start date
 *  - case_start_tod: case start time of day in seconds
 *  - input_file: config file (null terminated)
 *  - log_file: emulator log file (null terminated)
 *  - calendar: calendar identifier string
 */
struct EmulatorCreateConfig {
  int f_comm;
  int comp_id;
  int run_type;
  int start_ymd;
  int start_tod;
  int case_start_ymd;
  int case_start_tod;
  const char* input_file;
  const char* log_file;
  const char* calendar;
};

inline constexpr std::size_t EMULATOR_DIM_NAME_LENGTH = 64;

struct EmulatorDimensionDesc {
  char name[EMULATOR_DIM_NAME_LENGTH];
  int local_extent;
  int global_extent;
};

struct EmulatorGridDesc {
  int grid_type;
  int nx;
  int ny;

  int num_dimensions;
  const EmulatorDimensionDesc* dimensions;

  const int* col_gids;
  const double* lat;
  const double* lon;
  const double* area;
};


/**
 * @brief Fortran/C struct to create RegisteredField entry
 * Fields:
 * - role
 * - component
 * - attributes
 * - size
 * - data
 **/
struct RegisteredFieldDesc {
  const char* role;
  const char* component;
  size_t size;
};

} // extern "C"

#endif
