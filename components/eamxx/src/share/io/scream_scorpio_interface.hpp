#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "ekat/util/ekat_string_utils.hpp"
#include "share/field/field_tag.hpp"
#include "share/scream_types.hpp"

#include <vector>

/* C++/F90 bridge to F90 SCORPIO routines */

namespace scream {
namespace scorpio {

  using offset_t = std::int64_t;

  // WARNING: these values must match the ones of file_purpose_in and file_purpose_out
  // in the scream_scorpio_interface F90 module
  enum FileMode {
    Read = 1,
    Write = 2
  };
  /* All scorpio usage requires that the pio_subsystem is initialized. Happens only once per simulation */
  void eam_init_pio_subsystem(const int mpicom, const int atm_id = 0);
  /* Cleanup scorpio with pio_finalize */
  void eam_pio_finalize();
  /* Close a file currently open in scorpio */
  void eam_pio_closefile(const std::string& filename);
  /* Register a new file to be used for input/output with the scorpio module */
  void register_file(const std::string& filename, const FileMode mode);
  /* Sets the IO decompostion for all variables in a particular filename.  Required after all variables have been registered.  Called once per file. */
  void set_decomp(const std::string& filename);
  /* Sets the degrees-of-freedom for a particular variable in a particular file.  Called once for each variable, for each file. */
  void set_dof(const std::string &filename, const std::string &varname, const Int dof_len, const offset_t* x_dof);
  /* Register a dimension coordinate with a file. Called during the file setup. */
  void register_dimension(const std::string& filename,const std::string& shortname, const std::string& longname, const int length, const bool partitioned);
  /* Register a variable with a file.  Called during the file setup, for an output stream. */
  void register_variable(const std::string& filename, const std::string& shortname, const std::string& longname,
                         const std::string& units, const std::vector<std::string>& var_dimensions,
                         const std::string& dtype, const std::string& nc_dtype, const std::string& pio_decomp_tag);
  void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const std::string& meta_val);
  /* Register a variable with a file.  Called during the file setup, for an input stream. */
  void get_variable(const std::string& filename,const std::string& shortname, const std::string& longname,
                    const std::vector<std::string>& var_dimensions,
                    const std::string& dtype, const std::string& pio_decomp_tag);
  /* End the definition phase for a scorpio file.  Last thing called after all dimensions, variables, dof's and decomps have been set.  Called once per file.
   * Mandatory before writing or reading can happend on file. */
  void eam_pio_enddef(const std::string &filename);
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

extern "C" {
  /* Query whether the pio subsystem is inited or not */
  bool is_eam_pio_subsystem_inited();
  /* Checks if a file is already open, with the given mode */
  bool is_file_open_c2f(const char*&& filename, const int& mode);
  int get_int_attribute_c2f (const char*&& filename, const char*&& attr_name);
  void set_int_attribute_c2f (const char*&& filename, const char*&& attr_name, const int& value);
  int get_dimlen_c2f(const char*&& filename, const char*&& dimname);
  bool has_variable_c2f (const char*&& filename, const char*&& varname);
  /* Query a netCDF file for the time variable */
  double read_time_at_index_c2f(const char*&& filename, const int& time_index);
  double read_curr_time_c2f(const char*&& filename);
} // extern "C"

// The strings returned by e2str(const FieldTag&) are different from
// what existing nc files are already using. Besides upper/lower case
// differences, the column dimension (COL) is 'ncol' in nc files,
// but we'd like to keep 'COL' when printing our layouts, so we
// create this other mini helper function to get the name of a tag
// that is compatible with nc files. Note that tags that make no
// sense for an nc file are omitted. Namely, all those that have a
// field-dependent extent, such as vector dimensions. Those have to
// be "unpacked", storing a separate variable for each slice.

inline std::string get_nc_tag_name (const FieldTag& t, const int extent) {
  using namespace ShortFieldTagsNames;

  std::string name = "";
  switch(t) {
    case EL:
      name = "elem";
      break;
    case LEV:
      name = "lev";
      break;
    case ILEV:
      name = "ilev";
      break;
    case TL:
      name = "tl";
      break;
    case COL:
      name = "ncol";
      break;
    case GP:
      name = "gp";
      break;
    case CMP:
      name = "dim" + std::to_string(extent);
      break;
    // Added for rrtmgp - TODO revisit this paradigm, see comment in field_tag.hpp
    case NGAS:
      name = "ngas";
      break;
    case SWBND:
      name = "swband";
      break;
    case LWBND:
      name = "lwband";
      break;
    case SWGPT:
      name = "swgpt";
      break;
    case LWGPT:
      name = "lwgpt";
      break;
    default:
      EKAT_ERROR_MSG("Error! Field tag not supported in netcdf files.");
  }

  return name;
}

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
