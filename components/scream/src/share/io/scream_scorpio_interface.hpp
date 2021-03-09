#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "share/field/field_tag.hpp"
#include "share/scream_types.hpp"
#include <vector>

/* C++/F90 bridge to F90 SCORPIO routines */

// TODO, figure out a better way to define netCDF output type for fields
#ifdef SCREAM_CONFIG_IS_CMAKE
#  ifdef SCREAM_DOUBLE_PRECISION
  static constexpr int PIO_REAL = 6;
#  else
  static constexpr int PIO_REAL = 5;
#  endif // SCREAM_DOUBLE_PRECISION
#else // SCREAM_CONFIG_IS_CMAKE
  static constexpr int PIO_REAL = 6;
#endif // SCREAM_CONFIG_IS_CMAKE
static constexpr int PIO_INT = 4;

namespace scream {
namespace scorpio {
  /* All scorpio usage requires that the pio_subsystem is initialized. Happens only once per simulation */
  void eam_init_pio_subsystem(const int mpicom);
  /* Cleanup scorpio with pio_finalize */
  void eam_pio_finalize();
  /* Close a file currently open in scorpio */
  void eam_pio_closefile(const std::string& filename);
  /* Register a new file for output with scorpio module */
  void register_outfile(const std::string& filename);
  /* Register a new file to be used for input with the scorpio module */
  void register_infile(const std::string& filename);
  /* Every timestep each output file needs to be synced, call once per timestep, per file */
  void sync_outfile(const std::string& filename);
  /* Sets the IO decompostion for all variables in a particular filename.  Required after all variables have been registered.  Called once per file. */
  void set_decomp(const std::string& filename);
  /* Sets the degrees-of-freedom for a particular variable in a particular file.  Called once for each variable, for each file. */
  void set_dof(const std::string &filename, const std::string &varname, const Int dof_len, const Int* x_dof);
  /* Register a dimension coordinate with a file. Called during the file setup. */
  void register_dimension(const std::string& filename,const std::string& shortname, const std::string& longname, const int length);
  /* Register a variable with a file.  Called during the file setup, for an output stream. */
  void register_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  void register_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const std::vector<std::string>& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  /* Register a variable with a file.  Called during the file setup, for an input stream. */
  void get_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  void get_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const std::vector<std::string>& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  /* End the definition phase for a scorpio file.  Last thing called after all dimensions, variables, dof's and decomps have been set.  Called once per file.
   * Mandatory before writing or reading can happend on file. */
  void eam_pio_enddef(const std::string &filename);
  /* Called each timestep to update the timesnap for the last written output. */
  void pio_update_time(const std::string &filename, const Real time);

  /* Read data for a specific variable from a specific file. */
  void grid_read_data_array (const std::string &filename, const std::string &varname, const Int& dim_length, Real* hbuf);
  void grid_read_data_array (const std::string &filename, const std::string &varname, const std::vector<int>& dims, const Int& dim_length, const Int& padding, Real* hbuf);
  /* Write data for a specific variable to a specific file. */
  void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::vector<int>& dims, const Int& dim_length, const Int& padding, const Real* hbuf);

  /* Helper functions */
  void add_remove_padding(const int slow_dim_len, const int pad_dim_len, const int padding, const Real *hbuf_in, Real *hbuf_out, bool add_padding);
  void count_pio_atm_file();

extern "C" {
  /* Query whether the pio subsystem is inited or not */
  bool is_eam_pio_subsystem_inited();
  int  eam_pio_subsystem_comm ();
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

inline std::string get_nc_tag_name (const FieldTag& t) {
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
    default:
      EKAT_ERROR_MSG("Error! Field tag not supported in netcdf files.");
  }

  return name;
}

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
