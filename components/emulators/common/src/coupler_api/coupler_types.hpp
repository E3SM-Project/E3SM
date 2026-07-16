#ifndef COUPLER_TYPES
#define COUPLER_TYPES

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
  int  case_start_ymd;
  int  case_start_tod;
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
struct CouplingDesc {
  double* import_data;
  double* export_data;
  int     num_imports;
  int     num_exports;
  int     field_size;
};


struct CouplingStringList {
  const char* const* entries;
  int num_entries;
};

/**
  * @brief Description of Field designated for the coupler
  * Fields:
  * - longname: long name for field
  * - stdname: standardized name for field
  * - attrname: name for lookup in attribute vector
  * - units: physical units of field
  * - source: component that computes field
  * - destination: component that receives field
*/
struct CouplingFieldDesc {
  const char* longname;
  const char* stdname;
  const char* attrname;
  const char* units;
  const char* merge_type;
  CouplingStringList sources;
  CouplingStringList destination;
};

} //extern "C"

#endif
