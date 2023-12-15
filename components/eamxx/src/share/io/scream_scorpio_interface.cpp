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
  bool is_enddef_c2f(const char*&& filename);
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
void eam_init_pio_subsystem(const ekat::Comm& comm) {
  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  eam_init_pio_subsystem(fcomm);
}

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
int get_dimlen(const std::string& filename, const std::string& dimname)
{
  int ncid, dimid, err;
  PIO_Offset len;

  bool was_open = is_file_open_c2f(filename.c_str(),-1);
  if (not was_open) {
    register_file(filename,Read);
  }

  ncid = get_file_ncid_c2f (filename.c_str());
  err = PIOc_inq_dimid(ncid,dimname.c_str(),&dimid);
  EKAT_REQUIRE_MSG (err!=PIO_EBADDIM,
      "Error! Could not find dimension in the file.\n"
      " - filename : " + filename + "\n"
      " - dimname  : " + dimname + "\n"
      " - pio error: " + std::to_string(err) + "\n");
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving dimension id.\n"
      " - filename : " + filename + "\n"
      " - dimname  : " + dimname + "\n"
      " - pio error: " + std::to_string(err) + "\n");

  err = PIOc_inq_dimlen(ncid,dimid,&len);
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while querying dimension length.\n"
      " - filename : " + filename + "\n"
      " - dimname  : " + dimname + "\n"
      " - pio error: " + std::to_string(err) + "\n");

  if (not was_open) {
    eam_pio_closefile(filename);
  }

  return len;
}
/* ----------------------------------------------------------------- */
bool has_variable (const std::string& filename, const std::string& varname)
{
  int ncid, varid, err;

  bool was_open = is_file_open_c2f(filename.c_str(),-1);
  if (not was_open) {
    register_file(filename,Read);
  }

  ncid = get_file_ncid_c2f (filename.c_str());
  err = PIOc_inq_varid(ncid,varname.c_str(),&varid);
  if (err==PIO_ENOTVAR) {
    return false;
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving dimension id.\n"
      " - filename : " + filename + "\n"
      " - varname  : " + varname + "\n"
      " - pio error: " + std::to_string(err) + "\n");
  if (not was_open) {
    eam_pio_closefile(filename);
  }

  return true;
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
ekat::any get_any_attribute (const std::string& filename, const std::string& att_name) {
  register_file(filename,Read);
  auto ncid = get_file_ncid_c2f (filename.c_str());
  EKAT_REQUIRE_MSG (ncid>=0,
      "[get_any_attribute] Error! Could not retrieve file ncid.\n"
        " - filename : " + filename + "\n");

  int varid = PIO_GLOBAL;
  int err;

  nc_type type;
  PIO_Offset len;
  err = PIOc_inq_att(ncid,varid,att_name.c_str(),&type,&len);
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring global attribute.\n"
        " - filename : " + filename + "\n"
        " - attribute: " + att_name + "\n"
        " - pio error: " << err << "\n");

  EKAT_REQUIRE_MSG (len==1 || type==PIO_CHAR,
      "[get_any_attribute] Error! Only single value attributes allowed.\n"
        " - filename : " + filename + "\n"
        " - attribute: " + att_name + "\n"
        " - nc type  : " << type << "\n"
        " - att len  : " << len << "\n");

  ekat::any att;
  if (type==PIO_INT) {
    int val;
    err = PIOc_get_att(ncid,varid,att_name.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_DOUBLE) {
    double val;
    err = PIOc_get_att(ncid,varid,att_name.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_FLOAT) {
    float val;
    err = PIOc_get_att(ncid,varid,att_name.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_CHAR) {
    std::string val(len,'\0');
    err = PIOc_get_att(ncid,varid,att_name.c_str(),&val[0]);
    att.reset(val);
  } else {
    EKAT_ERROR_MSG ("[get_any_attribute] Error! Unsupported/unrecognized nc type.\n"
        " - filename : " + filename + "\n"
        " - attribute: " + att_name + "\n"
        " - nc type  : " << type << "\n");
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring global attribute.\n"
        " - filename : " + filename + "\n"
        " - attribute: " + att_name + "\n"
        " - pio error: " << err << "\n");

  eam_pio_closefile(filename);
  return att;
}
void set_any_attribute (const std::string& filename, const std::string& att_name, const ekat::any& att) {
  auto ncid = get_file_ncid_c2f (filename.c_str());
  int err;

  EKAT_REQUIRE_MSG (ncid>=0,
      "[set_any_attribute] Error! Could not retrieve file ncid.\n"
        " - filename : " + filename + "\n");

  bool redef = is_enddef_c2f(filename.c_str());
  if (redef) {
    err = PIOc_redef(ncid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "[set_any_attribute] Error! Something went wrong while re-opening def phase.\n"
          " - filename : " + filename + "\n"
          " - attribute: " + att_name + "\n"
          " - pio error: " << err << "\n");
  }

  int varid = PIO_GLOBAL;
  if (att.isType<int>()) {
    const int& data = ekat::any_cast<int>(att);
    err = PIOc_put_att(ncid,varid,att_name.c_str(),PIO_INT,1,&data);
  } else if (att.isType<double>()) {
    const double& data = ekat::any_cast<double>(att);
    err = PIOc_put_att(ncid,varid,att_name.c_str(),PIO_DOUBLE,1,&data);
  } else if (att.isType<float>()) {
    const float& data = ekat::any_cast<float>(att);
    err = PIOc_put_att(ncid,varid,att_name.c_str(),PIO_FLOAT,1,&data);
  } else if (att.isType<std::string>()) {
    const std::string& data = ekat::any_cast<std::string>(att);
    err = PIOc_put_att(ncid,varid,att_name.c_str(),PIO_CHAR,data.size(),data.data());
  } else {
    EKAT_ERROR_MSG ("[set_any_attribute] Error! Unsupported/unrecognized att type.\n"
        " - filename : " + filename + "\n"
        " - att name : " + att_name + "\n"
        " - att value: " << att << "\n"
        " - type info: " << att.content().type().name() << "\n");
  }

  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[set_any_attribute] Error! Something went wrong while setting global attribute.\n"
        " - filename : " + filename + "\n"
        " - attribute: " + att_name + "\n"
        " - pio error: " << err << "\n");

  if (redef) {
    err = PIOc_enddef(ncid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "[set_any_attribute] Error! Something went wrong while re-closing def phase.\n"
          " - filename : " + filename + "\n"
          " - attribute: " + att_name + "\n"
          " - pio error: " << err << "\n");
  }
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
void write_timestamp (const std::string& filename, const std::string& ts_name, const util::TimeStamp& ts)
{
  set_attribute(filename,ts_name,ts.to_string());
}
/* ----------------------------------------------------------------- */
util::TimeStamp read_timestamp (const std::string& filename, const std::string& ts_name)
{
  return util::str_to_time_stamp(get_attribute<std::string>(filename,ts_name));
}
/* ----------------------------------------------------------------- */
} // namespace scorpio
} // namespace scream
