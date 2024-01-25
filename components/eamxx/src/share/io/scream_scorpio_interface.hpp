#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "share/field/field_tag.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <vector>

/* C++/F90 bridge to F90 SCORPIO routines */

namespace scream {
namespace scorpio {

  using offset_t = std::int64_t;

  // WARNING: these values must match the ones of file_purpose_in and file_purpose_out
  // in the scream_scorpio_interface F90 module
  enum FileMode {
    Read = 1,
    Append = 2,
    Write = 4
  };
  /* All scorpio usage requires that the pio_subsystem is initialized. Happens only once per simulation */
  void eam_init_pio_subsystem(const ekat::Comm& comm);
  void eam_init_pio_subsystem(const int mpicom, const int atm_id = 0);
  /* Cleanup scorpio with pio_finalize */
  void eam_pio_finalize();
  /* Close a file currently open in scorpio */
  void eam_pio_closefile(const std::string& filename);
  void eam_flush_file(const std::string& filename);
  /* Register a new file to be used for input/output with the scorpio module */
  void register_file(const std::string& filename, const FileMode mode, int iotype);
  /* Sets the IO decompostion for all variables in a particular filename.  Required after all variables have been registered.  Called once per file. */
  int get_dimlen(const std::string& filename, const std::string& dimname);
  bool has_dim(const std::string& filename, const std::string& dimname);
  bool has_variable (const std::string& filename, const std::string& varname);
  bool has_attribute (const std::string& filename, const std::string& attname);
  bool has_attribute (const std::string& filename, const std::string& varname, const std::string& attname);
  void set_decomp(const std::string& filename);
  /* Sets the degrees-of-freedom for a particular variable in a particular file.  Called once for each variable, for each file. */
  void set_dof(const std::string &filename, const std::string &varname, const Int dof_len, const offset_t* x_dof);
  /* Register a dimension coordinate with a file. Called during the file setup. */
  void register_dimension(const std::string& filename,const std::string& shortname, const std::string& longname, const int length, const bool partitioned);
  /* Register a variable with a file.  Called during the file setup, for an output stream. */
  void register_variable(const std::string& filename, const std::string& shortname, const std::string& longname,
                         const std::string& units, const std::vector<std::string>& var_dimensions,
                         const std::string& dtype, const std::string& nc_dtype, const std::string& pio_decomp_tag);
  void register_variable(const std::string& filename, const std::string& shortname, const std::string& longname,
                         const std::vector<std::string>& var_dimensions,
                         const std::string& dtype, const std::string& pio_decomp_tag);
  void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const std::string& meta_val);
  void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const float meta_val);
  void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const double meta_val);
  void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, float& meta_val);
  void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, double& meta_val);
  void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, std::string& meta_val);
  /* Register a variable with a file.  Called during the file setup, for an input stream. */
  ekat::any get_any_attribute (const std::string& filename, const std::string& att_name);
  ekat::any get_any_attribute (const std::string& filename, const std::string& var_name, const std::string& att_name);
  void set_any_attribute (const std::string& filename, const std::string& att_name, const ekat::any& att);
  /* End the definition phase for a scorpio file.  Last thing called after all dimensions, variables, dof's and decomps have been set.  Called once per file.
   * Mandatory before writing or reading can happend on file. */
  void eam_pio_enddef(const std::string &filename);
  void eam_pio_redef(const std::string &filename);
  /* Called each timestep to update the timesnap for the last written output. */
  void pio_update_time(const std::string &filename, const double time);

  // Read data for a specific variable from a specific file. To read data that
  // isn't associated with a time index, or to read data at the most recent
  // time, set time_index to -1. Otherwise use the proper zero-based time index.
  template<typename T>
  void grid_read_data_array (const std::string &filename, const std::string &varname,
                             const int time_index, T* hbuf, const int buf_size);
  /* Write data for a specific variable to a specific file. */
  template<typename T>
  void grid_write_data_array(const std::string &filename, const std::string &varname,
                             const T* hbuf, const int buf_size);

  template<typename T>
  T get_attribute (const std::string& filename, const std::string& att_name)
  {
    auto att = get_any_attribute(filename,att_name);
    return ekat::any_cast<T>(att);
  }

  template<typename T>
  void set_attribute (const std::string& filename, const std::string& att_name, const T& att)
  {
    ekat::any a(att);
    set_any_attribute(filename,att_name,a);
  }

  // Shortcut to write/read to/from YYYYMMDD/HHMMSS attributes in the NC file
  void write_timestamp (const std::string& filename, const std::string& ts_name,
                        const util::TimeStamp& ts, const bool write_nsteps = false);
  util::TimeStamp read_timestamp (const std::string& filename,
                                  const std::string& ts_name,
                                  const bool read_nsteps = false);

extern "C" {
  /* Query whether the pio subsystem is inited or not */
  bool is_eam_pio_subsystem_inited();
  /* Checks if a file is already open, with the given mode */
  int get_file_ncid_c2f(const char*&& filename);
  // If mode<0, then simply checks if file is open, regardless of mode
  bool is_file_open_c2f(const char*&& filename, const int& mode);
  /* Query a netCDF file for the time variable */
  bool is_enddef_c2f(const char*&& filename);
  double read_time_at_index_c2f(const char*&& filename, const int& time_index);
  double read_curr_time_c2f(const char*&& filename);
  /* Query a netCDF file for the metadata associated w/ a variable */
  float get_variable_metadata_float_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name);
  double get_variable_metadata_double_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name);
} // extern "C"

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP
