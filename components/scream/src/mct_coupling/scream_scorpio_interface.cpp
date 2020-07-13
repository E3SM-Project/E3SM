#include "scream_scorpio_interface.hpp"
#include "scream_config.h"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include <string>

using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
  void register_infile_c(const char*&& filename);
  void set_decomp_c(const char*&& filename);
  void set_dof_c(const char*&& filename,const char*&& varname,const Int dof_len,const Int *x_dof);
  void grid_read_data_array_c_real(const char*&& filename, const char*&& varname, const Int dim1_length, Real *hbuf);
  void grid_read_data_array_c_int(const char*&& filename, const char*&& varname, const Int dim1_length, Int *hbuf);

  void grid_write_data_array_c_real_1d(const char*&& filename, const char*&& varname, const Int dim1_length, const Real* hbuf);
  void grid_write_data_array_c_real_2d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Real* hbuf);
  void grid_write_data_array_c_real_3d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Real* hbuf);
  void grid_write_data_array_c_real_4d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Real* hbuf);
  void grid_write_data_array_c_int_1d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int* hbuf);
  void grid_write_data_array_c_int_2d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int* hbuf);
  void grid_write_data_array_c_int_3d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int* hbuf);
  void grid_write_data_array_c_int_4d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Int* hbuf);
  void eam_init_pio_subsystem_c(const int mpicom, const int compid, const bool local);
  void eam_pio_finalize_c();
  void register_outfile_c(const char*&& filename);
  void sync_outfile_c(const char*&& filename);
  void eam_pio_closefile_c(const char*&& filename);
  void pio_update_time_c(const char*&& filename,const Real time);
  void register_dimension_c(const char*&& filename, const char*&& shortname, const char*&& longname, const int length);
  void register_variable_c(const char*&& filename,const char*&& shortname, const char*&& longname, const int numdims, const char** var_dimensions, const int dtype, const char*&& pio_decomp_tag);
  void get_variable_c(const char*&& filename,const char*&& shortname, const char*&& longname, const int numdims, const char** var_dimensions, const int dtype, const char*&& pio_decomp_tag);
  void eam_pio_enddef_c(const char*&& filename);

} // extern C

namespace scream {
namespace scorpio {
/* ----------------------------------------------------------------- */
void eam_init_pio_subsystem(const int mpicom, const int compid, const bool local) {
  eam_init_pio_subsystem_c(mpicom,compid,local);
}
/* ----------------------------------------------------------------- */
void eam_pio_finalize() {
  eam_pio_finalize_c();
}
/* ----------------------------------------------------------------- */
void register_outfile(const std::string& filename) {

  register_outfile_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void eam_pio_closefile(const std::string& filename) {

  eam_pio_closefile_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void register_infile(const std::string& filename) {

  register_infile_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void sync_outfile(const std::string& filename) {

  sync_outfile_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_decomp(const std::string& filename) {

  set_decomp_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_dof(const std::string& filename, const std::string& varname, const Int dof_len, const Int* x_dof) {

  set_dof_c(filename.c_str(),varname.c_str(),dof_len,x_dof);
}
/* ----------------------------------------------------------------- */
void pio_update_time(const std::string& filename, const Real time) {

  pio_update_time_c(filename.c_str(),time);
}
/* ----------------------------------------------------------------- */
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length) {

  register_dimension_c(filename.c_str(), shortname.c_str(), longname.c_str(), length);
}
/* ----------------------------------------------------------------- */
void get_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  get_variable_c(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions, dtype, pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  register_variable_c(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions, dtype, pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void eam_pio_enddef(const std::string &filename) {
  eam_pio_enddef_c(filename.c_str());
}
/* ----------------------------------------------------------------- */
void grid_read_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, Int *hbuf) {

  grid_read_data_array_c_int(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_read_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, Real *hbuf) {

  grid_read_data_array_c_real(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_1d(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_1d(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_1d(filename.c_str(),varname.c_str(),dim_length[0],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_2d(filename.c_str(),varname.c_str(),dim_length[0],dim_length[1],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_3d(filename.c_str(),varname.c_str(),dim_length[0],dim_length[1],dim_length[2],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_4d(filename.c_str(),varname.c_str(),dim_length[0],dim_length[1],dim_length[2],dim_length[3],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,1>& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_1d(filename.c_str(),varname.c_str(),dim_length[0],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_2d(filename.c_str(),varname.c_str(),dim_length[0],dim_length[1],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_3d(filename.c_str(),varname.c_str(),dim_length[0],dim_length[1],dim_length[2],hbuf);

};
/* ----------------------------------------------------------------- */

} // namespace scorpio
} // namespace scream
