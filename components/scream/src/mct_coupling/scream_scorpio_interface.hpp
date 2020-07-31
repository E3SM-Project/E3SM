#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

// TODO, figure out a better way to define netCDF output type for fields
#ifdef SCREAM_CONFIG_IS_CMAKE
#  ifdef SCREAM_DOUBLE_PRECISION
  const int PIO_REAL = 6;
#  else
  const int PIO_REAL = 5;
#  endif // SCREAM_DOUBLE_PRECISION
#else // SCREAM_CONFIG_IS_CMAKE
  const int PIO_REAL = 6;
#endif // SCREAM_CONFIG_IS_CMAKE
const int PIO_INT = 4;

using scream::Real;
using scream::Int;
namespace scream {
namespace scorpio {

  void eam_init_pio_subsystem(const int mpicom, const int compid, const bool local);
  void eam_pio_finalize();
  void eam_pio_closefile(const std::string& filename);
  void register_outfile(const std::string& filename);
  void register_infile(const std::string& filename);
  void sync_outfile(const std::string& filename);
  void set_decomp(const std::string& filename);
  void set_dof(const std::string &filename, const std::string &varname, const Int dof_len, const Int* x_dof);
  void register_dimension(const std::string& filename,const std::string& shortname, const std::string& longname, const int length);
  void register_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  void get_variable(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  void eam_pio_enddef(const std::string &filename);
  void pio_update_time(const std::string &filename, const Real time);

  void grid_read_data_array (const std::string &filename, const std::string &varname, const Int& dim_length, Real* hbuf);
  void grid_read_data_array (const std::string &filename, const std::string &varname, const Int& dim_length, Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Int* hbuf);

  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Real* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Int* hbuf);
  void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Int* hbuf);

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
