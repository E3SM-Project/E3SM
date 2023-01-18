#include "scream_scorpio_interface.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "scream_config.h"

#include "ekat/ekat_assert.hpp"
#include "share/scream_types.hpp"

#include <pio.h>

#include <string>


using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
  void register_file_c2f(const char*&& filename, const int& mode);
  void set_decomp_c2f(const char*&& filename);
  void set_dof_c2f(const char*&& filename,const char*&& varname,const Int dof_len,const std::int64_t *x_dof);
  void grid_read_data_array_c2f_int(const char*&& filename, const char*&& varname, const Int time_index, int *buf, const int buf_size);
  void grid_read_data_array_c2f_float(const char*&& filename, const char*&& varname, const Int time_index, float *buf, const int buf_size);
  void grid_read_data_array_c2f_double(const char*&& filename, const char*&& varname, const Int time_index, double *buf, const int buf_size);

  void grid_write_data_array_c2f_int(const char*&& filename, const char*&& varname, const int* buf, const int buf_size);
  void grid_write_data_array_c2f_float(const char*&& filename, const char*&& varname, const float* buf, const int buf_size);
  void grid_write_data_array_c2f_double(const char*&& filename, const char*&& varname, const double* buf, const int buf_size);
  void eam_init_pio_subsystem_c2f(const int mpicom, const int atm_id);
  void eam_pio_finalize_c2f();
  void eam_pio_closefile_c2f(const char*&& filename);
  void pio_update_time_c2f(const char*&& filename,const double time);
  void register_dimension_c2f(const char*&& filename, const char*&& shortname, const char*&& longname, const int global_length, const bool partitioned);
  void register_variable_c2f(const char*&& filename, const char*&& shortname, const char*&& longname,
                             const char*&& units, const int numdims, const char** var_dimensions,
                             const int dtype, const int nc_dtype, const char*&& pio_decomp_tag);
  void set_variable_metadata_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name, const char*&& meta_val);
  void get_variable_c2f(const char*&& filename,const char*&& shortname, const char*&& longname,
                        const int numdims, const char** var_dimensions,
                        const int dtype, const char*&& pio_decomp_tag);
  void eam_pio_enddef_c2f(const char*&& filename);
} // extern C

namespace scream {
namespace scorpio {

// Retrieve the int codes PIO uses to specify data types
int nctype (const std::string& type) {
  if (type=="int") {
    return PIO_INT;
  } else if (type=="float" || type=="single") {
    return PIO_FLOAT;
  } else if (type=="double") {
    return PIO_DOUBLE;
  } else if (type=="real") {
#if defined(SCREAM_DOUBLE_PRECISION)
    return PIO_DOUBLE;
#else
    return PIO_FLOAT;
#endif
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type '" + type + "'.\n");
  }
}
/* ----------------------------------------------------------------- */

void eam_init_pio_subsystem(const int mpicom, const int atm_id) {
  // TODO: Right now the compid has been hardcoded to 0 and the flag
  // to create a init a subsystem in SCREAM is hardcoded to true.
  // When surface coupling is established we will need to refactor this
  // routine to pass the appropriate values depending on if we are running
  // the full model or a unit test.
  eam_init_pio_subsystem_c2f(mpicom,atm_id);
}
/* ----------------------------------------------------------------- */
void eam_pio_finalize() {
  eam_pio_finalize_c2f();
}
/* ----------------------------------------------------------------- */
void register_file(const std::string& filename, const FileMode mode) {
  register_file_c2f(filename.c_str(),mode);
}
/* ----------------------------------------------------------------- */
void eam_pio_closefile(const std::string& filename) {

  eam_pio_closefile_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_decomp(const std::string& filename) {

  set_decomp_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_dof(const std::string& filename, const std::string& varname, const Int dof_len, const std::int64_t* x_dof) {

  set_dof_c2f(filename.c_str(),varname.c_str(),dof_len,x_dof);
}
/* ----------------------------------------------------------------- */
void pio_update_time(const std::string& filename, const double time) {

  pio_update_time_c2f(filename.c_str(),time);
}
/* ----------------------------------------------------------------- */
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length, const bool partitioned) {

  register_dimension_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), length, partitioned);
}
/* ----------------------------------------------------------------- */
void get_variable(const std::string &filename, const std::string& shortname, const std::string& longname,
                  const std::vector<std::string>& var_dimensions,
                  const std::string& dtype, const std::string& pio_decomp_tag) {

  /* Convert the vector of strings that contains the variable dimensions to a char array */
  const int numdims = var_dimensions.size();
  std::vector<const char*> var_dimensions_c(numdims);
  for (int ii = 0;ii<numdims;++ii) 
  {
    var_dimensions_c[ii] = var_dimensions[ii].c_str();
  }
  get_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(),
                   numdims, var_dimensions_c.data(), nctype(dtype), pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname,
                       const std::string& units, const std::vector<std::string>& var_dimensions,
                       const std::string& dtype, const std::string& nc_dtype, const std::string& pio_decomp_tag) {

  /* Convert the vector of strings that contains the variable dimensions to a char array */
  const int numdims = var_dimensions.size();
  std::vector<const char*> var_dimensions_c(numdims);
  for (int ii = 0;ii<numdims;++ii) 
  {
    var_dimensions_c[ii] = var_dimensions[ii].c_str();
  }
  register_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(),
                        units.c_str(), numdims, var_dimensions_c.data(),
                        nctype(dtype), nctype(nc_dtype), pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const std::string& meta_val) {
  set_variable_metadata_c2f(filename.c_str(),varname.c_str(),meta_name.c_str(),meta_val.c_str());
}
/* ----------------------------------------------------------------- */
void eam_pio_enddef(const std::string &filename) {
  eam_pio_enddef_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
template<>
void grid_read_data_array<int>(const std::string &filename, const std::string &varname,
                          const int time_index, int *hbuf, const int buf_size) {
  grid_read_data_array_c2f_int(filename.c_str(),varname.c_str(),time_index,hbuf,buf_size);
}
template<>
void grid_read_data_array<float>(const std::string &filename, const std::string &varname,
                                const int time_index, float *hbuf, const int buf_size) {
  grid_read_data_array_c2f_float(filename.c_str(),varname.c_str(),time_index,hbuf,buf_size);
}
template<>
void grid_read_data_array<double>(const std::string &filename, const std::string &varname,
                                  const int time_index, double *hbuf, const int buf_size) {
  grid_read_data_array_c2f_double(filename.c_str(),varname.c_str(),time_index,hbuf,buf_size);
}
/* ----------------------------------------------------------------- */
template<>
void grid_write_data_array<int>(const std::string &filename, const std::string &varname, const int* hbuf, const int buf_size) {
  grid_write_data_array_c2f_int(filename.c_str(),varname.c_str(),hbuf,buf_size);
}
template<>
void grid_write_data_array<float>(const std::string &filename, const std::string &varname, const float* hbuf, const int buf_size) {
  grid_write_data_array_c2f_float(filename.c_str(),varname.c_str(),hbuf,buf_size);
}
template<>
void grid_write_data_array<double>(const std::string &filename, const std::string &varname, const double* hbuf, const int buf_size) {
  grid_write_data_array_c2f_double(filename.c_str(),varname.c_str(),hbuf,buf_size);
}
/* ----------------------------------------------------------------- */
} // namespace scorpio
} // namespace scream
