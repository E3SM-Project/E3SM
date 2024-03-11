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
  void register_file_c2f(const char*&& filename, const int& mode, const int& iotype);
  int get_file_mode_c2f(const char*&& filename);
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
  void eam_pio_flush_file_c2f(const char*&& filename);
  void pio_update_time_c2f(const char*&& filename,const double time);
  void register_dimension_c2f(const char*&& filename, const char*&& shortname, const char*&& longname, const int global_length, const bool partitioned);
  void register_variable_c2f(const char*&& filename, const char*&& shortname, const char*&& longname,
                             const char*&& units, const int numdims, const char** var_dimensions,
                             const int dtype, const int nc_dtype, const char*&& pio_decomp_tag);
  void set_variable_metadata_char_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name, const char*&& meta_val);
  void set_variable_metadata_float_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name, const float meta_val);
  void set_variable_metadata_double_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name, const double meta_val);
  float get_variable_metadata_float_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name);
  double get_variable_metadata_double_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name);
  void get_variable_metadata_char_c2f (const char*&& filename, const char*&& varname, const char*&& meta_name, char*&& meta_val);
  void eam_pio_redef_c2f(const char*&& filename);
  void eam_pio_enddef_c2f(const char*&& filename);
  bool is_enddef_c2f(const char*&& filename);
} // extern C

namespace scream {
namespace scorpio {

// Retrieve the int codes PIO uses to specify data types
int nctype (const std::string& type) {
  if (type=="int") {
    return PIO_INT;
  } else if (type=="int64") {
    return PIO_INT64;
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
std::string nctype2str (const int type) {
  for (auto t : {"int", "int64", "float", "double"}) {
    if (nctype(t)==type) return t;
  }
  return "UNKNOWN";
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
void register_file(const std::string& filename, const FileMode mode, const int iotype) {
  register_file_c2f(filename.c_str(),mode,iotype);
}
/* ----------------------------------------------------------------- */
void eam_pio_closefile(const std::string& filename) {

  eam_pio_closefile_c2f(filename.c_str());
}
void eam_flush_file(const std::string& filename) {
  eam_pio_flush_file_c2f(filename.c_str());
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
bool has_dim (const std::string& filename, const std::string& dimname)
{
  int ncid, dimid, err;

  bool was_open = is_file_open_c2f(filename.c_str(),-1);
  if (not was_open) {
    register_file(filename,Read);
  }

  ncid = get_file_ncid_c2f (filename.c_str());
  err = PIOc_inq_dimid(ncid,dimname.c_str(),&dimid);
  if (err==PIO_EBADDIM) {
    return false;
  }

  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving dimension id.\n"
      " - filename : " + filename + "\n"
      " - dimname  : " + dimname + "\n"
      " - pio error: " + std::to_string(err) + "\n");
  if (not was_open) {
    eam_pio_closefile(filename);
  }

  return true;
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
      "Error! Something went wrong while retrieving variable id.\n"
      " - filename : " + filename + "\n"
      " - varname  : " + varname + "\n"
      " - pio error: " + std::to_string(err) + "\n");
  if (not was_open) {
    eam_pio_closefile(filename);
  }

  return true;
}

bool has_attribute (const std::string& filename, const std::string& attname)
{
  return has_attribute(filename,"GLOBAL",attname);
}

bool has_attribute (const std::string& filename, const std::string& varname, const std::string& attname)
{
  int ncid, varid, attid, err;

  bool was_open = is_file_open_c2f(filename.c_str(),-1);
  if (not was_open) {
    register_file(filename,Read);
  }

  // Get file id
  ncid = get_file_ncid_c2f (filename.c_str());

  // Get var id
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    err = PIOc_inq_varid(ncid,varname.c_str(),&varid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong while retrieving variable id.\n"
        " - filename : " + filename + "\n"
        " - varname  : " + varname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
  }

  // Get att id
  err = PIOc_inq_attid(ncid,varid,attname.c_str(),&attid);
  if (err==PIO_ENOTATT) {
    return false;
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving attribute id.\n"
      " - filename : " + filename + "\n"
      " - varname  : " + varname + "\n"
      " - attname  : " + attname + "\n"
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
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length, const bool partitioned)
{
  int mode = get_file_mode_c2f(filename.c_str());
  std::string mode_str = mode==Read ? "Read" : (mode==Write ? "Write" : "Append");
  if (mode!=Write) {
    // Ensure the dimension already exists, and that it has the correct size (if not unlimited)
    int ncid,dimid,unlimid,err;

    ncid = get_file_ncid_c2f (filename.c_str());
    err = PIOc_inq_dimid(ncid,shortname.c_str(),&dimid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Could not retrieve dimension id from file open in " + mode_str + " mode.\n"
        " - filename: " + filename + "\n"
        " - dimension : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");

    err = PIOc_inq_unlimdim(ncid,&unlimid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong querying for the unlimited dimension id.\n"
        " - filename: " + filename + "\n"
        " - dimension : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
    if (length==0) {
      EKAT_REQUIRE_MSG ( unlimid==dimid,
          "Error! Input dimension is unlimited, but does not appear to be unlimited in the file (open in " + mode_str + " mode).\n"
          " - filename: " + filename + "\n"
          " - dimension : " + shortname + "\n"
          " - pio error: " + std::to_string(err) + "\n");
    } else {
      EKAT_REQUIRE_MSG ( unlimid!=dimid,
          "Error! Input dimension is not unlimited, but it appears to be unlimited in the file (open in " + mode_str + " mode).\n"
          " - filename: " + filename + "\n"
          " - dimension : " + shortname + "\n"
          " - pio error: " + std::to_string(err) + "\n");

      int len_from_file = get_dimlen(filename,shortname);
      EKAT_REQUIRE_MSG (length==len_from_file,
          "Error! Input dimension length does not match the one from the file  (open in " + mode_str + " mode).\n"
          " - filename: " + filename + "\n"
          " - dimension : " + shortname + "\n"
          " - input dim length: " + std::to_string(length) + "\n"
          " - file dim length : " + std::to_string(len_from_file) + "\n");
    }
  }

  register_dimension_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), length, partitioned);
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string& filename, const std::string& shortname, const std::string& longname,
                       const std::vector<std::string>& var_dimensions,
                       const std::string& dtype, const std::string& pio_decomp_tag)
{
  // This overload does not require to specify an nc data type, so it *MUST* be used when the
  // file access mode is either Read or Append. Either way, a) the var should be on file already,
  // and b) so should be the dimensions
  EKAT_REQUIRE_MSG (has_variable(filename,shortname),
      "Error! This overload of register_variable *assumes* the variable is already in the file, but wasn't found.\n"
      " - filename: " + filename + "\n"
      " - varname : " + shortname + "\n");
  for (const auto& dimname : var_dimensions) {
    int len = get_dimlen (filename,dimname);
    // WARNING! If the dimension was not yet registered, it will be registered as a NOT partitioned dim.
    //          If this dim should be partitioned, then register it *before* the variable
    register_dimension(filename,dimname,dimname,len,false);
  }
  register_variable(filename,shortname,longname,"",var_dimensions,dtype,"",pio_decomp_tag);
}
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname,
                       const std::string& units_in, const std::vector<std::string>& var_dimensions,
                       const std::string& dtype, const std::string& nc_dtype_in, const std::string& pio_decomp_tag)
{
  // Local copies, since we can modify them in case of defaults
  auto units = units_in;
  auto nc_dtype = nc_dtype_in;

  int mode = get_file_mode_c2f(filename.c_str());
  std::string mode_str = mode==Read ? "Read" : (mode==Write ? "Write" : "Append");

  bool has_var = has_variable(filename,shortname);
  if (mode==Write) {
    EKAT_REQUIRE_MSG ( units!="" and nc_dtype!="",
        "Error! Missing valid units and/or nc_dtype arguments for file open in Write mode.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n");
  } else {
    EKAT_REQUIRE_MSG ( has_var,
        "Error! Variable not found in file open in " + mode_str + " mode.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n");
  }


  if (has_var) {
    // The file already exists or the var was already registered.
    // Make sure we're registering the var with the same specs
    int ncid = get_file_ncid_c2f (filename.c_str());
    int vid,ndims,err;
    err = PIOc_inq_varid(ncid,shortname.c_str(),&vid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong retrieving variable id.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
    err = PIOc_inq_varndims(ncid,vid,&ndims);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong inquiring the number of dimensions of a variable.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
    std::vector<int> dims(ndims);
    err = PIOc_inq_vardimid(ncid,vid,dims.data());
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong inquiring the dimensions ids of a variable.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
    std::vector<std::string> dims_from_file(ndims);
    for (int i=0; i<ndims; ++i) {
      dims_from_file[i].resize(PIO_MAX_NAME);
      err =  PIOc_inq_dimname(ncid,dims[i],&dims_from_file[i][0]);
      // IMPORTANT: must remove trailing null chars, to get *correct* string size
      dims_from_file[i] = ekat::trim(dims_from_file[i],'\0');
      EKAT_REQUIRE_MSG (err==PIO_NOERR,
          "Error! Something went wrong inquiring dimension name.\n"
          " - filename: " + filename + "\n"
          " - dimid   : " + std::to_string(dims[i]) + "\n"
          " - pio error: " + std::to_string(err) + "\n");
    }

    // Here, let's only try to access var_dimensions[0] when we know for sure
    // that var_dimensions is actually dimensioned (i.e., .size()>0)
    if (var_dimensions.size()>0) {
      if (mode==Read && (dims_from_file[0]=="time" && var_dimensions[0]!="time")) {
        // For Read operations, we may not consider "time" as a field dimension, so if the
        // input file has "time", simply disregard it in this check.
        dims_from_file.erase(dims_from_file.begin());
      }
    } else {
      if (mode==Read && (dims_from_file[0]=="time")) {
        // For Read operations, we may not consider "time" as a field dimension, so if the
        // input file has "time", simply disregard it in this check.
        dims_from_file.erase(dims_from_file.begin());
      }      
    }
    std::reverse(dims_from_file.begin(),dims_from_file.end());
    EKAT_REQUIRE_MSG(var_dimensions==dims_from_file,
        "Error! Input variable dimensions do not match the ones from the file.\n"
        " - filename  : " + filename + "\n"
        " - varname   : " + shortname + "\n"
        " - input dims: (" + ekat::join(var_dimensions,",") + ")\n"
        " - file dims : (" + ekat::join(dims_from_file,",") + ")\n");

    // Check nc dtype only if user bothered specifying it (if mode=Read, probably the user doesn't care)
    nc_type type;
    err = PIOc_inq_vartype(ncid,vid,&type);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong while inquiring variable data type.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n"
        " - pio error: " + std::to_string(err) + "\n");
    EKAT_REQUIRE_MSG (nc_dtype=="" or type==nctype(nc_dtype),
        "Error! Input NC data type does not match the one from the file.\n"
        " - filename: " + filename + "\n"
        " - varname : " + shortname + "\n"
        " - input dtype : " + nc_dtype + "\n"
        " - file dtype  : " + nctype2str(type) + "\n");
    nc_dtype = nctype2str(type);

    // Get var units (if set)
    int natts;
    err = PIOc_inq_varnatts(ncid,vid,&natts);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong while inquiring a variable's number of attributes.\n"
          " - filename: " + filename + "\n"
          " - varname : " + shortname + "\n"
          " - pio error: " + std::to_string(err) + "\n");
    std::string units_from_file(PIO_MAX_NAME,'\0');
    std::string att_name(PIO_MAX_NAME,'\0');
    for (int i=0; i<natts; ++i) {
      err = PIOc_inq_attname(ncid,vid,i,&att_name[0]);
      EKAT_REQUIRE_MSG (err==PIO_NOERR,
          "Error! Something went wrong while inquiring the name of a variable attribute.\n"
            " - filename: " + filename + "\n"
            " - varname : " + shortname + "\n"
            " - pio error: " + std::to_string(err) + "\n");
      if (att_name=="units") {
        err = PIOc_get_att_text(ncid,vid,"units",&units_from_file[0]);
        EKAT_REQUIRE_MSG (err==PIO_NOERR,
            "Error! Something went wrong while retrieving a variable units.\n"
              " - filename: " + filename + "\n"
              " - varname : " + shortname + "\n"
              " - pio error: " + std::to_string(err) + "\n");
        EKAT_REQUIRE_MSG (units=="" or units==units_from_file,
            "Error! Input variable units do not match the ones from the file.\n"
            " - filename: " + filename + "\n"
            " - varname : " + shortname + "\n"
            " - input units: " + units + "\n"
            " - file units : " + units_from_file + "\n");

        // If mode=Read and user didn't bother with units, but units _are_ in the file, then
        // make sure we call register_variable with correct units (to avoid later checks failures)
        units = units_from_file;
      }
    }
  }

  // Convert the vector of strings that contains the variable dimensions to a char array
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
void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const float meta_val) {
  set_variable_metadata_float_c2f(filename.c_str(),varname.c_str(),meta_name.c_str(),meta_val);
}
/* ----------------------------------------------------------------- */
void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const double meta_val) {
  set_variable_metadata_double_c2f(filename.c_str(),varname.c_str(),meta_name.c_str(),meta_val);
}
/* ----------------------------------------------------------------- */
void set_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, const std::string& meta_val) {
  set_variable_metadata_char_c2f(filename.c_str(),varname.c_str(),meta_name.c_str(),meta_val.c_str());
}
/* ----------------------------------------------------------------- */
void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, float& meta_val) {
  meta_val = get_variable_metadata_float_c2f(filename.c_str(),varname.c_str(),meta_name.c_str());
}
/* ----------------------------------------------------------------- */
void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, double& meta_val) {
  meta_val = get_variable_metadata_double_c2f(filename.c_str(),varname.c_str(),meta_name.c_str());
}
/* ----------------------------------------------------------------- */
void get_variable_metadata (const std::string& filename, const std::string& varname, const std::string& meta_name, std::string& meta_val) {
  meta_val.resize(256);
  get_variable_metadata_char_c2f(filename.c_str(),varname.c_str(),meta_name.c_str(),&meta_val[0]);

  // If terminating char is not found, meta_val simply uses all 256 chars
  if (meta_val.find('\0')!=std::string::npos) {
    meta_val.resize(meta_val.find('\0'));
  }
}
/* ----------------------------------------------------------------- */
ekat::any get_any_attribute (const std::string& filename, const std::string& att_name) {
  auto out = get_any_attribute(filename,"GLOBAL",att_name);
  return out;
}
/* ----------------------------------------------------------------- */
ekat::any get_any_attribute (const std::string& filename, const std::string& var_name, const std::string& att_name) {
  register_file(filename,Read);
  auto ncid = get_file_ncid_c2f (filename.c_str());
  EKAT_REQUIRE_MSG (ncid>=0,
      "[get_any_attribute] Error! Could not retrieve file ncid.\n"
        " - filename : " + filename + "\n");

  int varid;
  int err;
  if (var_name=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    err = PIOc_inq_varid(ncid, var_name.c_str(), &varid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "[get_any_attribute] Error! Something went wrong while inquiring variable id.\n"
          " - filename : " + filename + "\n"
          " - variable : " + var_name + "\n"
          " - attribute: " + att_name + "\n"
          " - pio error: " << err << "\n");
  }

  nc_type type;
  PIO_Offset len;
  err = PIOc_inq_att(ncid,varid,att_name.c_str(),&type,&len);
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring global attribute.\n"
        " - filename : " + filename + "\n"
        " - variable : " + var_name + "\n"
        " - attribute: " + att_name + "\n"
        " - pio error: " << err << "\n");

  EKAT_REQUIRE_MSG (len==1 || type==PIO_CHAR,
      "[get_any_attribute] Error! Only single value attributes allowed.\n"
        " - filename : " + filename + "\n"
        " - variable : " + var_name + "\n"
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
        " - variable : " + var_name + "\n"
        " - attribute: " + att_name + "\n"
        " - nc type  : " << type << "\n");
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring global attribute.\n"
        " - filename : " + filename + "\n"
        " - variable : " + var_name + "\n"
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
void eam_pio_redef(const std::string &filename) {
  eam_pio_redef_c2f(filename.c_str());
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
void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps)
{
  set_attribute(filename,ts_name,ts.to_string());
  if (write_nsteps) {
    set_attribute(filename,ts_name+"_nsteps",ts.get_num_steps());
  }
}
/* ----------------------------------------------------------------- */
util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps)
{
  auto ts = util::str_to_time_stamp(get_attribute<std::string>(filename,ts_name));
  if (read_nsteps and has_attribute(filename,ts_name+"_nsteps")) {
    ts.set_num_steps(get_attribute<int>(filename,ts_name+"_nsteps"));
  }
  return ts;
}
/* ----------------------------------------------------------------- */
} // namespace scorpio
} // namespace scream
