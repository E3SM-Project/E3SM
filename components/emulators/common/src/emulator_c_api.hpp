/**
 * @file emulator_c_api.hpp
 * @brief Define structs and functions to hide all implementation details from Fortran API via opaque pointers
*/
#ifndef EMULATOR_C_API
#define EMULATOR_C_API
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
*  - input_file: config file (null terminated)
*  - log_file: emulator log file (null terminated)
*/
struct EmulatorCreateConfig {
  int  f_comm;
  int  comp_id;
  int  run_type;
  int  start_ymd;
  int  start_tod;
  const char* input_file;
  const char* log_file;
};

/**
 * @brief Description for the grid decomposition
 * 
 * Fields:
* - grid_type: structured/unstructured
* - nx: 
* - ny
* - num_local_cols
* - num_global_cols
* - col_gids
* - lat
* - lon
* - area
*/
struct EmulatorGridDesc {
  int grid_type;
  int nx;
  int ny;
  int num_local_cols;
  int num_global_cols;
  const int*    col_gids;
  const double* lat;
  const double* lon;
  const double* area;
};

/**
 * @brief Description of import and export fields to/from the coupler
 * Fields:
 *  - import_data
 *  - export_data
 *  - num_imports
 *  - num_exports
 *  - field_size
*
*/
struct EmulatorCouplingDesc {
  double* import_data;
  double* export_data;
  int     num_imports;
  int     num_exports;
  int     field_size;
};

/// Opaque handle type in C/Fortran:
/// actually points to an EmulatorComp in C++.
void* emulator_create(const char* kind,
                      const EmulatorCreateConfig* cfg);

void  emulator_set_grid_data(void* handle,
                             const EmulatorGridDesc* grid);

void  emulator_setup_coupling(void* handle,
                              EmulatorCouplingDesc* cpl);

void  emulator_init(void* handle);
void  emulator_run(void* handle, int dt);
void  emulator_finalize(void* handle);
void  emulator_print_info(void* handle);

} // extern "C"
#endif
