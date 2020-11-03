#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

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
  void grid_read_data_array (const std::string &filename, const std::string &varname, const Int& dim_length, Int* hbuf);
  /* Write data for a specific variable to a specific file. */
  void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Int* hbuf);
  /* Read and write data routines that allow for passage of multi dimensional arrays.  May no longer be needed, but keep for now just in case. */
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Int* hbuf);

  void count_pio_atm_file();
} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
