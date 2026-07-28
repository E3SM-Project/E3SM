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
  int  f_comm;
  int  comp_id;
  int  run_type;
  int  start_ymd;
  int  start_tod;
  int  case_start_ymd;
  int  case_start_tod;
  const char* input_file;
  const char* log_file;
  const char* calendar;
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
  * @brief Metadata of Field that may be coupled
  * Fields:
  * - longname: long name for field
  * - stdname: standardized name for field
  * - attrname: name for lookup in attribute vector
*/
struct FieldAttributes{
  const char* name;
  const char* long_name;
  const char* standard_name;
  const char* units;
} ;

/**
 * @brief Fortran/C struct to create RegisteredField entry
 * Fields:
 * - component
 * - attributes
 * - size
 * - data
 **/
struct RegisteredFieldDesc {
  const char* role;
  const char* component;
  FieldAttributes attributes;
  size_t size;
  double* data;
};

/**
  * @brief Description of How Field is Coupled
  * Fields:
  * - merged_type: method to merge data from multiple sources
  * - source: component that computes field
  * - destination: component that receives field
*/
struct CouplingFieldDesc {
  const char* merge_type;
  CouplingStringList sources;
  CouplingStringList destinations;
};

} //extern "C"

#endif
