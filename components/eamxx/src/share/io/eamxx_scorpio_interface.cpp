#include "eamxx_scorpio_interface.hpp"
#include "eamxx_shr_interface_c2f.hpp"

#include "eamxx_config.h"

#include <ekat/ekat_assert.hpp>
#include <ekat/util/ekat_string_utils.hpp>

#include <pio.h>

#include <numeric>

namespace scream {
namespace scorpio {

// This class is an implementation detail, and therefore it is hidden inside
// a cpp file. All customers of IO capabilities must use the common interfaces
// exposed in the header file of this source file.
// This class simply serves as a container for persistent IO data.
struct ScorpioSession
{
public:
  static ScorpioSession& instance () {
    static ScorpioSession s;
    return s;
  }

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  strmap_t<PIOFile>                     files;
  strmap_t<std::shared_ptr<PIODecomp>>  decomps;

  // In the above map, we would like to label decomps as dtype_dim1$N1_dim2$N2...
  // where N$i is the global length of dim$i. However, it *may* happen that we use
  // two different decompositions for the same global layout, which would clash
  // the names. For this reason, we append at the end an increasing counter,
  // which disambiguate between globally-equivalent partitions.
  // When adding a new decomp, we check this map to see if another decomp
  // already exists with the same global layout. If so, we first check to see
  // if that decomp is equivalent to the new one *on all ranks*. If yes, we
  // recycle it, otherwise we create a new PIO decomp.
  // strmap_t<std::list<std::string>>    decomp_global_layout_to_decomp_name;
  // strmap_t<int>                       decomp_global_layout_to_counter;

  int         pio_sysid        = -1;
  int         pio_type_default = -1;
  int         pio_rearranger   = -1;
  int         pio_format       = -1;

  ekat::Comm  comm;

private:

  ScorpioSession () = default;
};

// --------------------------------------------------------------------------------------------- //

template<typename S, typename D>
void copy_data (const S* src, D* dst, int n) {
  for (int i=0; i<n; ++i) {
    dst[i] = src[i];
  }
}

// Utility for common IO operation failure
void check_scorpio_noerr (const int err,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: PIOc_" + pioc_func_name + "\n");
}

void check_scorpio_noerr (const int err, const std::string& filename,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - filename: " + filename + "\n" 
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: " + pioc_func_name + "\n");
}

void check_scorpio_noerr (const int err,
                          const std::string& filename,
                          const std::string& entity_type,
                          const std::string& entity_name,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - filename: " + filename + "\n" 
      " - " + entity_type + ": " + entity_name + "\n" 
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: " + pioc_func_name + "\n");
}

// Return name of a shared ptr to PIO entity (to use inside ekat::join)
std::string get_entity_name (const std::shared_ptr<const PIOEntity>& e)
{
  return e->name;
}

template<typename T>
std::string print_map_keys (const std::map<std::string,T>& map) {
  std::string s;
  for (const auto& it : map) {
    s += it.first + ",";
  }
  s.pop_back();
  return s;
}

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
  } else if (type=="char") {
    return PIO_CHAR;
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type '" + type + "'.\n");
  }
  return -1;
}

template<typename T>
int nctype () {
  if (std::is_same<T,int>::value) {
    return PIO_INT;
  } else if (std::is_same<T,std::int64_t>::value) {
    return PIO_INT64;
  } else if (std::is_same<T,float>::value) {
    return PIO_FLOAT;
  } else if (std::is_same<T,double>::value) {
    return PIO_DOUBLE;
  } else if (std::is_same<T,std::string>::value) {
    return PIO_CHAR;
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type.\n");
  }
}

template<typename T>
const void* ncdata (const T& v) {
  return reinterpret_cast<const void*>(&v);
}
template<>
const void* ncdata (const std::string& v) {
  return reinterpret_cast<const void*>(v.data());
}

template<typename T>
PIO_Offset nclen (const T& v) {
  return 1;
}
template<>
PIO_Offset nclen (const std::string& v) {
  return v.size();
}

std::string refine_dtype (const std::string& dtype) {
  if (dtype=="real") {
#if defined(SCREAM_DOUBLE_PRECISION)
    return "double";
#else
    return "float";
#endif
  } else if (dtype=="single") {
    return "float";
  } else {
    return dtype;
  }
}

size_t dtype_size (const std::string& dtype) {
  if (dtype=="int") {
    return sizeof(int);
  } else if (dtype=="int64") {
    return sizeof(long long);
  } else if (dtype=="float") {
    return sizeof(float);
  } else if (dtype=="double") {
    return sizeof(double);
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type '" + dtype + "'.\n");
  }
  return 0;
}

int pio_iotype (IOType iotype) {
  int iotype_int;
  auto& s = ScorpioSession::instance();
  switch(iotype){
    case IOType::DefaultIOType: iotype_int = s.pio_type_default;                    break;
    case IOType::NetCDF:        iotype_int = static_cast<int>(PIO_IOTYPE_NETCDF);   break;
    case IOType::PnetCDF:       iotype_int = static_cast<int>(PIO_IOTYPE_PNETCDF);  break;
    case IOType::Adios:         iotype_int = static_cast<int>(PIO_IOTYPE_ADIOS);    break;
    case IOType::Adiosc:        iotype_int = static_cast<int>(PIO_IOTYPE_ADIOSC);   break;
    case IOType::Hdf5:          iotype_int = static_cast<int>(PIO_IOTYPE_HDF5);     break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported iotype.\n");
  }
  return iotype_int;
}

// ====================== Local utilities ========================== //

namespace impl {
// Note: these utilities are used in this file to retrieve PIO entities,
//       so that we implement all checks once (rather than in every function)

// Small struct that allows to quickly open a file (in Read mode) if it wasn't open.
// If the file had to be open, when the struct is deleted, it will release the file.
struct PeekFile {
  PeekFile(const std::string& filename_in) {
    filename = filename_in;
    was_open = is_file_open(filename);
    if (not was_open) {
      register_file(filename,Read);
    }
    auto& s = ScorpioSession::instance();
    file = &s.files[filename];
  }

  ~PeekFile () {
    if (not was_open) {
      // Note: this function _could_ throw, but it should not happen (unless someone
      //       else called release_file twice). That's b/c we either are not the only
      //       customer (so nothing to be done other than a ref count decrement) or
      //       the file was open in Read mode, in which case it doesn't need to do
      //       much in scorpio.
      release_file(filename);
    }
  }

  const PIOFile*  file;
  std::string     filename;
  bool            was_open;
};

PIOFile& get_file (const std::string& filename,
                   const std::string& context)
{
  auto& s = ScorpioSession::instance();

  EKAT_REQUIRE_MSG (s.files.count(filename)==1,
      "Error! Could not retrieve the file. File not open.\n"
      " - filename: " + filename + "\n"
      "Context:\n"
      " " + context + "\n");

  return s.files.at(filename);
}

PIODim& get_dim (const std::string& filename,
                 const std::string& dimname,
                 const std::string& context)
{
  const auto& f = get_file(filename,context);
  EKAT_REQUIRE_MSG (f.dims.count(dimname)==1,
      "Error! Could not retrieve dimension. Dimension not found.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n"
      " - dims on file: " + print_map_keys(f.dims) + "\n"
      "Context:\n"
      " " + context + "\n");

  return *f.dims.at(dimname);
}

PIOVar& get_var (const std::string& filename,
                 const std::string& varname,
                 const std::string& context)
{
  const auto& f = get_file(filename,context);
  EKAT_REQUIRE_MSG (f.vars.count(varname)==1,
      "Error! Could not retrieve variable. Variable not found.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n"
      " - vars on file : " + print_map_keys(f.vars) + "\n"
      "Context:\n"
      " " + context + "\n");

  return *f.vars.at(varname);
}

} // namespace impl

// ====================== Global IO operations ======================= // 

void init_subsystem(const ekat::Comm& comm, const int atm_id)
{
  auto& s = ScorpioSession::instance();
  s.comm = comm;

  EKAT_REQUIRE_MSG (s.pio_sysid==-1,
      "Error! Attmept to re-initialize pio subsystem.\n");

#ifdef SCREAM_CIME_BUILD
  s.pio_sysid        = shr_get_iosysid_c2f(atm_id);
  s.pio_type_default = shr_get_iotype_c2f(atm_id);
  s.pio_rearranger   = shr_get_rearranger_c2f(atm_id);
  s.pio_format       = shr_get_ioformat_c2f(atm_id);
#else
  // Use some reasonable defaults for standalone EAMxx tests
  int stride = 1;
  int base = 0;

  s.pio_rearranger = PIO_REARR_SUBSET;
  s.pio_format     = PIO_64BIT_DATA;
#if PIO_USE_PNETCDF
  s.pio_type_default = pio_iotype(IOType::PnetCDF);
#elif PIO_USE_NETCDF
  s.pio_type_default = pio_iotype(IOType::NetCDF);
#else
#error "Standalone EAMxx requires either PNETCDF or NETCDF iotype to be available in Scorpio"
#endif

  auto err = PIOc_Init_Intracomm(comm.mpi_comm(), comm.size(), stride, base, s.pio_rearranger, &s.pio_sysid);
  check_scorpio_noerr (err,"init_subsystem", "Init_Intracomm");

  // Unused in standalone mode
  (void) atm_id;
#endif

  static_assert (sizeof(offset_t)==sizeof(PIO_Offset),
      "Error! PIO was configured with PIO_OFFSET not a 64-bit int.\n");
}

bool is_subsystem_inited () {
  return ScorpioSession::instance().pio_sysid!=-1;
}

void finalize_subsystem ()
{
  auto& s = ScorpioSession::instance();

  // TODO: should we simply return instead? I think trying to finalize twice
  //       *may* be a sign of possible bugs, though with Catch2 testing
  //       I *think* there may be some issue with how the code is run.
  EKAT_REQUIRE_MSG (s.pio_sysid!=-1,
      "Error! PIO subsystem was already finalized.\n");

  for (auto& it : s.files) {
    EKAT_REQUIRE_MSG (it.second.num_customers==0,
      "Error! ScorpioSession::finalize called, but a file is still in use elsewhere.\n"
      " - filename: " + it.first + "\n");
  }
  s.files.clear();

  for (auto& it : s.decomps) {
    EKAT_REQUIRE_MSG (it.second.use_count()==1,
      "Error! ScorpioSession::finalize called, but a decomp is still stored elsewhere.\n"
      " - decomp name: " + it.first + "\n");

    int err = PIOc_freedecomp(s.pio_sysid,it.second->ncid);
    check_scorpio_noerr(err,"finalize_subsystem","freedecomp");
  }
  s.decomps.clear();

#ifndef SCREAM_CIME_BUILD
  // Don't finalize in CIME builds, since the coupler will take care of it
  PIOc_finalize (s.pio_sysid);
#endif

  s.pio_sysid        = -1;
  s.pio_type_default = -1;
  s.pio_format       = -1;
  s.pio_rearranger   = -1;
}

// ========================= File operations ===================== //

void register_file (const std::string& filename,
                    const FileMode mode,
                    const IOType iotype)
{
  auto& s = ScorpioSession::instance();
  auto& f = s.files[filename];
  EKAT_REQUIRE_MSG (f.mode==Unset || f.mode==mode,
      "Error! File was already opened with a different mode.\n"
      " - filename: " + filename + "\n"
      " - old mode: " + e2str(f.mode) + "\n"
      " - new mode: " + e2str(mode) + "\n");
  EKAT_REQUIRE_MSG (f.mode==Unset || f.iotype==iotype,
      "Error! File was already opened with a different iotype.\n"
      " - filename: " + filename + "\n"
      " - old type: " + iotype2str(f.iotype) + "\n"
      " - new type: " + iotype2str(iotype) + "\n");

  if (f.mode == Unset) {
    // First time we ask for this file. Call PIO open routine(s)
    int err;
    int iotype_int = pio_iotype(iotype);
    if (mode & Read) {
      auto write = mode & Write ? PIO_WRITE : PIO_NOWRITE;
      err = PIOc_openfile(s.pio_sysid,&f.ncid,&iotype_int,filename.c_str(),write);
      f.enddef = true;
    } else {
      err = PIOc_createfile(s.pio_sysid,&f.ncid,&iotype_int,filename.c_str(),s.pio_format);
      f.enddef = false;
    }

    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong while opening a file.\n"
        " - filename : " + filename + "\n"
        " - file mode: " + e2str(mode) + "\n"
        " - pio error: " + std::to_string(err) + "\n");

    f.mode = mode;
    f.iotype = iotype;
    f.name = filename;

    if (mode & Read) {
      // Read all dims/vars from file
      PIO_Offset len;
      int ndims, nvars, ngatts, unlimdimid;
      err = PIOc_inq(f.ncid, &ndims, &nvars, &ngatts, &unlimdimid);
      check_scorpio_noerr(err,f.name,"register_file","inq");

      char name[PIO_MAX_NAME];
      for (int idim=0; idim<ndims; ++idim) {
        err = PIOc_inq_dim(f.ncid,idim,name,&len);
        check_scorpio_noerr (err,f.name,"register_file","inq_dim");
        auto dim = f.dims[name] = std::make_shared<PIODim>();
        dim->name = name;
        dim->fid = f.ncid;
        dim->length = len;
        dim->ncid = idim;
        dim->unlimited = idim==unlimdimid;

        if (dim->unlimited) {
          f.time_dim = dim;
        }
      }

      auto find_dim = [&](int dimid) -> std::shared_ptr<const PIODim> {
        std::shared_ptr<const PIODim> d;
        for (auto it : f.dims) {
          if (it.second->ncid==dimid) {
            d = it.second;
          }
        }
        EKAT_REQUIRE_MSG (d!=nullptr,
            "Error! Could not locat dimension id in the file.\n"
            " - filename: " + f.name + "\n"
            " - dim id  : " + std::to_string(dimid) + "\n");
        return d;
      };
      int dtype,natts;
      int dimids[PIO_MAX_DIMS];
      for (int ivar=0; ivar<nvars; ++ivar) {
        err = PIOc_inq_var(f.ncid,ivar,name,PIO_MAX_NAME,&dtype,&ndims,dimids,&natts);
        check_scorpio_noerr(err,f.name,"register_file","inq_var");

        auto var = f.vars[name] = std::make_shared<PIOVar>();
        var->name = name;
        var->ncid = ivar;
        var->fid = f.ncid;
        for (int idim=0; idim<ndims; ++idim) {
          auto dim = find_dim(dimids[idim]);
          if (dim->unlimited) {
            var->time_dep = true;
            var->num_records = dim->length;
          } else {
            var->dims.push_back(find_dim(dimids[idim]));
          }
        }

        for (const auto& dt : {"int","int64","float","double","char"}) {
          if (nctype(dt)==dtype) {
            var->dtype = var->nc_dtype = dt;
            break;
          }
        }
        EKAT_REQUIRE_MSG (var->dtype!="",
            "Error! Variable data type not supported.\n"
            " - filename: " + filename + "\n"
            " - varname : " + var->name + "\n"
            " - nc dtype: " + std::to_string(dtype) + "\n");

        for (int iatt=0; iatt<natts; ++iatt) {
          err = PIOc_inq_attname(f.ncid,ivar,iatt,name);
          check_scorpio_noerr(err,f.name,"register_file","inq_attname");
          err = PIOc_inq_attlen(f.ncid,ivar,name,&len);
          check_scorpio_noerr(err,f.name,"register_file","inq_attlen");
          if (std::string(name)=="units") {
            err = PIOc_get_att_text(f.ncid,ivar,"units",name);
            name[len] = '\0';
            check_scorpio_noerr(err,f.name,"register_file","get_att_text");
            var->units = name;
            break;
          }
        }
      }
    }
  }
  ++f.num_customers;
}

void release_file  (const std::string& filename)
{
  auto& f = impl::get_file(filename,"scorpio::release_file");

  --f.num_customers;
  if (f.num_customers>0) {
    return;
  }

  int err;
  if (f.mode & Write) {
    err = PIOc_sync(f.ncid);
    check_scorpio_noerr (err,f.name,"release_file","sync");
  }

  err = PIOc_closefile(f.ncid);
  check_scorpio_noerr (err,f.name,"release_file","closefile");

  auto& s = ScorpioSession::instance();
  s.files.erase(filename);
}

void flush_file (const std::string &filename)
{
  auto& f = impl::get_file(filename,"scorpio::sync_file");
  
  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Cannot call sync_file. File is read-only.\n"
      " - filename: " + filename + "\n");

  int err = PIOc_sync(f.ncid);
  check_scorpio_noerr (err,f.name,"sync_file","sync");
}

void redef(const std::string &filename)
{
  auto& f = impl::get_file(filename,"scorpio::redef");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not call redef on the input file. File is read-only.\n"
      " - filename: " + filename + "\n");

  if (f.enddef) {
    int err = PIOc_redef(f.ncid);
    check_scorpio_noerr (err,f.name,"redef","redef");
    f.enddef = false;
  }
}

void enddef(const std::string &filename)
{
  auto& f = impl::get_file(filename,"scorpio::enddef");

  if (not f.enddef) {
    int err = PIOc_enddef(f.ncid);
    check_scorpio_noerr (err,f.name,"enddef","enddef");
    f.enddef = true;
  }
}

bool is_file_open (const std::string& filename, const FileMode mode)
{
  auto& s = ScorpioSession::instance();
  auto it = s.files.find(filename);
  if (it==s.files.end()) return false;

  return mode==Unset || (mode & it->second.mode);
}

// =================== Dimensions operations ======================= //

void define_dim (const std::string& filename, const std::string& dimname, const int length)
{
  auto& f = impl::get_file(filename,"scorpio::define_dim");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not define dimension. File is read-only.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  auto& dim = f.dims[dimname];

  bool unlimited = length==0;

  if (dim==nullptr) {
    EKAT_REQUIRE_MSG (f.mode!=Append,
        "Error! Cannot add a new dim when the file is open in append mode.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n");
    // Create new dimension
    dim = std::make_shared<PIODim>();
    dim->name = dimname;
    dim->fid = f.ncid;
    dim->length = length;
    dim->unlimited = unlimited;

    // Define the dimension in PIO
    int err = PIOc_def_dim(f.ncid,dimname.c_str(),dim->length,&dim->ncid);
    check_scorpio_noerr (err,f.name,"dimension",dimname,"define_dim","def_dim");
  } else {
    // Already defined. Check that the dim specs are the same.
    EKAT_REQUIRE_MSG (unlimited==dim->unlimited,
        "Error! Redefining dimension with different unlimited flag.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - old unlimited:" + (dim->unlimited ? "yes" : "no" )+ "\n"
        " - new unlimited:" + (unlimited ? "yes" : "no") + "\n");

    EKAT_REQUIRE_MSG (unlimited || length==dim->length,
        "Error! Redefining dimension with a different (local) length.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - old length:" + std::to_string(dim->length)+ "\n"
        " - new length:" + std::to_string(length) + "\n");
  }
}

bool has_dim (const std::string& filename,
              const std::string& dimname,
              const int length)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  auto it = pf.file->dims.find(dimname);
  if (it==pf.file->dims.end()) {
    return false;
  } else if (length==0) {
    return it->second->unlimited;
  } else if (length>0) {
    return it->second->length==length;
  }
  return true;
}

int get_dimlen (const std::string& filename, const std::string& dimname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  EKAT_REQUIRE_MSG (has_dim(filename,dimname),
      "Error! Could not inquire dimension length. The dimension is not in the file.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  return pf.file->dims.at(dimname)->length;
}

int get_dimlen_local (const std::string& filename, const std::string& dimname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  EKAT_REQUIRE_MSG (has_dim(filename,dimname),
      "Error! Could not inquire dimension local length. The dimension is not in the file.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  const auto& dim = pf.file->dims.at(dimname);
  return dim->offsets==nullptr ? dim->length : dim->offsets->size();
}

bool has_time_dim (const std::string& filename)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  return pf.file->time_dim!=nullptr;
}

int get_time_len (const std::string& filename)
{
  EKAT_REQUIRE_MSG (has_time_dim(filename),
      "Error! Could not inquire time dimension length. The time dimension is not in the file.\n"
      " - filename: " + filename + "\n");

  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  return pf.file->time_dim->length;
}

std::string get_time_name (const std::string& filename)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  EKAT_REQUIRE_MSG (pf.file->time_dim!=nullptr,
      "Error! Could not inquire time dimension name. The time dimension is not in the file.\n"
      " - filename: " + filename + "\n");

  return pf.file->time_dim->name;
}

void reset_time_dim_len (const std::string& filename, const int new_length)
{
  EKAT_REQUIRE_MSG (has_time_dim(filename),
      "Error! Could not reset time dimension length. The time dimension is not in the file.\n"
      " - filename: " + filename + "\n");

  auto& f = impl::get_file(filename,"scorpio::reset_time_dim_len");

  // Reset dim length
  EKAT_REQUIRE_MSG (new_length<f.time_dim->length,
      "Error! New time dimension length must be shorter than the current one.\n"
      "  - file name: " + filename + "\n"
      "  - curr len : " + std::to_string(f.time_dim->length) + "\n"
      "  - new len  : " + std::to_string(new_length) + "\n");
  f.time_dim->length = new_length;

  // Reset number of records counter for each time dep var
  for (auto it : f.vars) {
    auto& v = *it.second;
    if (v.time_dep) {
      v.num_records = new_length;
    }
  }
}

// =================== Decompositions operations ==================== //

// NOTES:
//  - this is a local function, we don't expose it. It's only called inside other scorpio utilities
//  - we don't really *need* filename, it's only to print more context in case of errors
void set_var_decomp (PIOVar& var,
                     const std::string& filename)
{
  for (size_t i=1; i<var.dims.size(); ++i) {
    EKAT_REQUIRE_MSG (var.dims[i]->offsets==nullptr,
        "Error! We currently only allow decomposition on slowest-striding dimension.\n"
        "       Generalizing is not complicated, but it was not a priority.\n"
        " - filename: " + filename + "\n"
        " - varname : " + var.name + "\n"
        " - var dims: " + ekat::join(var.dims,get_entity_name,",") + "\n"
        " - bad dim : " + var.dims[i]->name + "\n");
  }
  EKAT_REQUIRE_MSG (var.dims[0]->offsets!=nullptr,
      "Error! Calling set_var_decomp, but the var first dimension does not appear to be decomposed.\n"
      " - filename: " + filename + "\n"
      " - varname : " + var.name  + "\n"
      " - var dims: " + ekat::join(var.dims,get_entity_name,",") + "\n");
  EKAT_REQUIRE_MSG (var.decomp==nullptr,
      "Error! You should have invalidated var.decomp before attempting to reset it.\n"
      " - filename  : " + filename + "\n"
      " - varname   : " + var.name  + "\n"
      " - var decomp: " + var.decomp->name  + "\n");

  // Create decomp name: dtype-dim1<len1>_dim2<len2>_..._dimk<lenN>
  std::shared_ptr<const PIODim> decomp_dim;
  std::string decomp_tag = var.dtype + "-";
  for (auto d : var.dims) {
    decomp_tag += d->name + "<" + std::to_string(d->length) + ">_";
  }
  decomp_tag.pop_back(); // remove trailing underscore

  // Check if a decomp with this name already exists
  auto& s = ScorpioSession::instance();
  auto& decomp = s.decomps[decomp_tag];
#ifndef NDEBUG
  // Extra check: all ranks must agree on whether they have the decomposition!
  // If they don't agree, some rank will be stuck in a PIO call, waiting for others
  int found = decomp==nullptr ? 0 : 1;
  int min_found, max_found;
  const auto& comm = ScorpioSession::instance().comm;
  comm.all_reduce(&found,&min_found,1,MPI_MIN);
  comm.all_reduce(&found,&max_found,1,MPI_MAX);
  EKAT_REQUIRE_MSG(min_found==max_found,
      "Error! Decomposition already present on some ranks but not all.\n"
      " - filename: " + filename + "\n"
      " - varname : " + var.name + "\n"
      " - var dims: " + ekat::join(var.dims,get_entity_name,",") + "\n"
      " - decopm tag: " + decomp_tag + "\n");
#endif

  if (decomp==nullptr) {
    // We haven't create this decomp yet. Go ahead and create one
    decomp = std::make_shared<PIODecomp>();
    decomp->name = decomp_tag;
    decomp->dim = var.dims[0];

    int ndims = var.dims.size();

    // Get ALL dims global lengths, and compute prod of *non-decomposed* dims
    std::vector<int> gdimlen = {decomp->dim->length};
    int non_decomp_dim_prod = 1;
    for (int idim=1; idim<ndims; ++idim) {
      auto d = var.dims[idim];
      gdimlen.push_back(d->length);
      non_decomp_dim_prod *= d->length;
    }

    // Create offsets list
    const auto& dim_offsets = *decomp->dim->offsets;
    int dim_loc_len = dim_offsets.size();
    decomp->offsets.resize (non_decomp_dim_prod*dim_loc_len);
    for (int idof=0; idof<dim_loc_len; ++idof) {
      auto dof_offset = dim_offsets[idof];
      auto beg = decomp->offsets.begin()+ idof*non_decomp_dim_prod;
      auto end = beg + non_decomp_dim_prod;
      std::iota (beg,end,non_decomp_dim_prod*dof_offset);
    }

    // Create PIO decomp
    int maplen = decomp->offsets.size();
    PIO_Offset* compmap = reinterpret_cast<PIO_Offset*>(decomp->offsets.data());
    int err = PIOc_init_decomp(s.pio_sysid,nctype(var.dtype),ndims,gdimlen.data(),
                               maplen,compmap, &decomp->ncid,s.pio_rearranger,
                               nullptr,nullptr);

    check_scorpio_noerr(err,filename,"decomp",decomp_tag,"set_var_decomp","InitDecomp");
  }

  // Set decomp data in the var
  var.decomp = decomp;
}

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const std::vector<offset_t>& my_offsets,
                     const bool allow_reset)
{
  auto& s = ScorpioSession::instance();
  auto& f = impl::get_file(filename,"scorpio::set_decomp");
  auto& dim = impl::get_dim(filename,dimname,"scorpio::set_dim_decomp");

  EKAT_REQUIRE_MSG (not dim.unlimited,
      "Error! Cannot partition an unlimited dimension.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  if (dim.offsets!=nullptr) {
    if (allow_reset) {
      // We likely won't need the previously created decomps that included this dimension.
      // So, as we remove decomps from vars that have this dim, keep track of their name,
      // so that we can free them later *if no other users of them remain*.
      std::set<std::string> decomps_to_remove;
      for (auto it : f.vars) {
        auto v = it.second;
        if (v->decomp!=nullptr and v->decomp->dim->name==dimname) {
          decomps_to_remove.insert(v->decomp->name);
          v->decomp = nullptr;
        }
      }
      for (const auto& dn : decomps_to_remove) {
        if (s.decomps.at(dn).use_count()==1) {
          auto decomp = s.decomps.at(dn);
          // There is no other customer of this decomposition, so we can safely free it
          int err = PIOc_freedecomp(s.pio_sysid,decomp->ncid);
          check_scorpio_noerr(err,filename,"decomp",dn,"set_dim_decomp","freedecomp");
          s.decomps.erase(dn);
        }
      }
    } else {
      // Check that the offsets are (globally) the same
      int same = *dim.offsets==my_offsets;
      const auto& comm = ScorpioSession::instance().comm;
      comm.all_reduce(&same,1,MPI_MIN);
      EKAT_REQUIRE_MSG(same==1,
          "Error! Attempt to redefine a decomposition with a different dofs distribution.\n"
          " - filename: " + filename + "\n"
          " - dimname : " + dimname + "\n"
          "If you are attempting to redefine the decomp, call this function with throw_if_changing_decomp=false.\n");

      // Same decomposition, so we can just return
      return;
    }
  }
  
  // Check that offsets are less than the global dimension length
  for (auto o : my_offsets) {
    EKAT_REQUIRE_MSG (o>=0 && o<dim.length,
        "Error! Offset for dimension decomposition is out of bounds.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - dim glen: " + std::to_string(dim.length) + "\n"
        " - offset  : " + std::to_string(o) + "\n");
  }

  dim.offsets = std::make_shared<std::vector<offset_t>>(my_offsets);

  // If vars were already defined, we need to process them,
  // and create the proper PIODecomp objects.
  for (auto it : f.vars) {
    if (ekat::contains(it.second->dim_names(),dimname)) {
      set_var_decomp (*it.second,filename);
    }
  }
}

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const offset_t start, const offset_t count,
                     const bool allow_reset)
{
  std::vector<offset_t> offsets(count);
  std::iota(offsets.begin(),offsets.end(),start);
  set_dim_decomp(filename,dimname,offsets,allow_reset);
}

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const bool allow_reset)
{
  const auto& comm = ScorpioSession::instance().comm;

  const int glen = get_dimlen(filename,dimname);
  int len = glen / comm.size();
  if (comm.rank() < (glen % comm.size())) {
    ++len;
  }

  offset_t offset = len;
  comm.scan(&offset,1,MPI_SUM);
  offset -= len; // scan is inclusive, but we need exclusive

  set_dim_decomp (filename,dimname,offset,len,allow_reset);
}

// ================== Variable operations ================== //

// Define var on output file (cannot call on Read/Append files)
void define_var (const std::string& filename, const std::string& varname,
                 const std::string& units, const std::vector<std::string>& dimensions,
                 const std::string& dtype, const std::string& nc_dtype,
                 const bool time_dep)
{
  auto& f = impl::get_file(filename,"scorpio::define_var");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not define variable. File is read-only.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  EKAT_REQUIRE_MSG (not time_dep || f.time_dim!=nullptr,
      "Error! Cannot define time-dependent variable: no time dimension defined.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  if (f.vars.count(varname)==0) {
    EKAT_REQUIRE_MSG (f.mode!=Append,
        "Error! Cannot add a new var when the file is open in append mode.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n");
    // Create new variable
    auto var = std::make_shared<PIOVar>();
    var->name = varname;
    var->fid = f.ncid;
    var->units = units;
    var->dtype = refine_dtype(dtype);
    var->nc_dtype = refine_dtype(nc_dtype);
    var->time_dep = time_dep;
    int ndims = dimensions.size() + (time_dep ? 1 : 0);
    std::vector<int> dimids;
    if (time_dep) {
      dimids.push_back(f.time_dim->ncid);
    }
    for (const auto& dname : dimensions) {
      EKAT_REQUIRE_MSG (has_dim(filename,dname),
          "Error! Cannot create variable. Dimension not found.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - var dims : " + ekat::join(dimensions,",") + "\n"
          " - file dims: " + print_map_keys(f.dims) + "\n");
      auto dim = f.dims.at(dname);
      var->dims.push_back(dim);
      dimids.push_back(dim->ncid);
    }

    // Define the variable in PIO
    int err = PIOc_def_var(f.ncid,varname.c_str(),nctype(nc_dtype),ndims,dimids.data(),&var->ncid);
    check_scorpio_noerr(err,f.name,"variable",varname,"define_var","def_var");

    f.vars[varname] = var;

    if (units!="") {
      // Add units attribute
      set_attribute(filename,varname,"units",units);
    }

    if (var->dims.size()>0 and var->dims[0]->offsets!=nullptr) {
      set_var_decomp (*var,filename);
    }
  } else {
    const auto& var = f.vars.at(varname);
    // The variable was already defined. Check that important metadata is the same
    EKAT_REQUIRE_MSG (var->units==units,
        "Error! Attempt to redefine variable with different units.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old units: " + var->units + "\n"
          " - new units: " +      units + "\n");
    EKAT_REQUIRE_MSG (var->dtype==refine_dtype(dtype),
        "Error! Attempt to redefine variable with different data type.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old dtype: " + var->dtype + "\n"
          " - new dtype: " +      refine_dtype(dtype) + "\n");
    EKAT_REQUIRE_MSG (var->nc_dtype==refine_dtype(nc_dtype),
        "Error! Attempt to redefine variable with different PIO data type.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old pio dtype: " + var->nc_dtype + "\n"
          " - new pio dtype: " +      refine_dtype(nc_dtype) + "\n");
    EKAT_REQUIRE_MSG (var->time_dep==time_dep,
        "Error! Attempt to redefine variable with different time dep flag.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old time_dep: " + (var->time_dep ? "yes" : "no") + "\n"
          " - new time_dep: " + (     time_dep ? "yes" : "no") + "\n");
    const auto var_dims = ekat::join(var->dims,get_entity_name,",");
    EKAT_REQUIRE_MSG (var_dims==ekat::join(dimensions,","),
        "Error! Attempt to redefine variable with different dimensions.\n"
          " - filename: " + filename + "\n"
          " - varname : " + varname + "\n"
          " - old dims: " + var_dims + "\n"
          " - new dims: " + ekat::join(dimensions,",") + "\n");
  }
}

void define_var (const std::string& filename, const std::string& varname,
                 const std::vector<std::string>& dimensions,
                 const std::string& dtype,
                 const bool time_dependent)
{
  define_var(filename,varname,"",dimensions,dtype,dtype,time_dependent);
}

// This overload is not exposed externally. Also, filename is only
// used to print it in case there are errors
void change_var_dtype (PIOVar& var,
                       const std::string& dtype,
                       const std::string& filename)
{
  if (refine_dtype(dtype)==refine_dtype(var.dtype)) {
    // The type is not changing, nothing to do
    return;
  }

  var.dtype = refine_dtype(dtype);
  if (var.decomp) {
    // Re-decompose the variable, with new data type
    var.decomp = nullptr;
    set_var_decomp (var,filename);
  }
}

void change_var_dtype (const std::string& filename,
                       const std::string& varname,
                       const std::string& dtype)
{
  auto& var = impl::get_var(filename,varname,"scorpio::change_var_dtype");
  change_var_dtype(var,dtype,filename);
}

bool has_var (const std::string& filename, const std::string& varname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  return pf.file->vars.count(varname)==1;
}

const PIOVar& get_var (const std::string& filename,
                       const std::string& varname)
{
  return impl::get_var(filename,varname,"scorpio::get_var");
}

void define_time (const std::string& filename, const std::string& units, const std::string& time_name)
{
  auto& f = impl::get_file(filename,"scorpio::define_time");
  EKAT_REQUIRE_MSG (f.time_dim==nullptr,
      "Error! Attempt to redeclare unlimited dimension.\n"
      " - filename: " + filename + "\n");

  define_dim(filename,time_name,0);
  f.time_dim = f.dims.at(time_name);

  define_var(filename,time_name,units,{},"double","double",true);
}

void mark_dim_as_time (const std::string& filename, const std::string& dimname)
{
  auto& f = impl::get_file(filename,"scorpio::mark_dim_as_time");

  EKAT_REQUIRE_MSG (not has_time_dim(filename),
      "Error! Resetting the time dimension is not allowed once set (even if it's the same).\n"
      " - filename: " + filename + "\n"
      " - old time dim name: " + f.time_dim->name + "\n"
      " - new time dim name: " + dimname + "\n");

  EKAT_REQUIRE_MSG (f.mode==Read,
      "Error! Cannot interpret dimension as 'time' dim. File not in Read mode.\n"
      " - filename : " + filename + "\n"
      " - file mode: " + e2str(f.mode) + "\n");

  if (f.time_dim==nullptr) {
    EKAT_REQUIRE_MSG (has_dim(filename,dimname),
        "Error! Cannot interpret dimension as 'time' dim. Dimension not found.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n");

    auto dim = f.dims.at(dimname);
    f.time_dim = dim;

    // If a var has "time" in its dims (must be the 1st dim!),
    // remove it. Recall that we only store non-time dims in
    // the list of var dims.
    for (auto& it : f.vars) {
      auto& v = it.second;
      if (v->dims.size()>0 and v->dims[0]->name==dimname) {
        v->dims.erase(v->dims.begin());
        v->size = -1;
        v->time_dep = true;
      }
    }
  } else {
    EKAT_REQUIRE_MSG (f.time_dim->name==dimname,
        "Error! Attempt to change the time dimension.\n"
        " - filenama    : " + filename + "\n"
        " - old time dim: " + f.time_dim->name + "\n"
        " - new time dim: " + dimname + "\n");
  }
}

// Update value of time variable, increasing time dim length
void update_time(const std::string &filename, const double time) {
  const auto& f = impl::get_file(filename,"scorpio::update_time");
        auto& time_dim = *f.time_dim;
  const auto& var = impl::get_var(filename,time_dim.name,"scorpio::update_time");

  PIO_Offset index = time_dim.length;
  int err = PIOc_put_var1(f.ncid,var.ncid,&index,&time);
  check_scorpio_noerr (err,f.name,"update time","put_var1");
  ++time_dim.length;
}

double get_time (const std::string& filename, const int time_index)
{
  impl::PeekFile pf (filename);

  const auto& time_name = pf.file->time_dim->name;
  double t;
  read_var(filename,time_name,&t,time_index);
  return t;
}

std::vector<double> get_all_times (const std::string& filename)
{
  impl::PeekFile pf (filename);
  const auto& dim = *pf.file->time_dim;

  std::vector<double> times (dim.length);
  for (int i=0; i<dim.length; ++i) {
    read_var (filename, dim.name, &times[i], i);
  }
  return times;
}

// Read variable into user provided buffer.
// If time dim is present, read given time slice (time_index=-1 means "read last record).
// If time dim is not present, time_index must be -1 (error out otherwise)
template<typename T>
void read_var (const std::string &filename, const std::string &varname, T* buf, const int time_index)
{
  EKAT_REQUIRE_MSG (buf!=nullptr,
      "Error! Cannot read from provided pointer. Invalid buffer pointer.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  const auto& f = impl::get_file(filename,"scorpio::read_var");
        auto& var = impl::get_var(filename,varname,"scorpio::read_var");

  // If the input pointer type already matches var.dtype, this is a no-op
  change_var_dtype(var,get_dtype<T>(),filename);

  int err;
  int frame = -1;
  if (var.time_dep) {
    frame = time_index>=0 ? time_index : f.time_dim->length-1;
    EKAT_REQUIRE_MSG (frame<f.time_dim->length,
        "Error! Time index out of bounds.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n"
        " - time idx: " + std::to_string(time_index) + "\n"
        " - time len: " + std::to_string(f.time_dim->length));
    err = PIOc_setframe(f.ncid,var.ncid,frame);
    check_scorpio_noerr (err,f.name,"variable",varname,"read_var","setframe");
  }

  std::string pioc_func;
  if (var.decomp) {
    // A decomposed variable, requires read_darray
    err = PIOc_read_darray(f.ncid,var.ncid,var.decomp->ncid,var.decomp->offsets.size(),buf);
    pioc_func = "read_darray";
  } else {
    // A non-decomposed variable, use PIOc_get_var(a)

    // If nc data type doesn't match the input pointer, we need to use the var internal buffer
    void* io_buf = buf;
    if (var.dtype!=var.nc_dtype) {
      if (var.size==-1) {
        var.size = 1;
        for (auto d : var.dims) {
          var.size *= d->length;
        }
        var.buf.resize(var.size*dtype_size(var.nc_dtype));
      }
      io_buf = var.buf.data();
    }

    if (frame>=0) {
      // We need to get the start/count for each dimension
      int ndims = var.dims.size();
      std::vector<PIO_Offset> start (ndims+1,0), count(ndims+1); // +1 for time
      start[0] = frame;
      count[0] = 1;
      for (int idim=0; idim<ndims; ++idim) {
        count[idim+1] = var.dims[idim]->length;
      }
      err = PIOc_get_vara(f.ncid,var.ncid,start.data(),count.data(),io_buf);
      pioc_func = "get_vara";
    } else {
      err = PIOc_get_var(f.ncid,var.ncid,io_buf);
      pioc_func = "get_var";
    }

    // If we used the var tmp buffer, copy back into the user-provided pointer
    if (var.dtype!=var.nc_dtype) {
      if (var.nc_dtype=="int") {
        copy_data(reinterpret_cast<int*>(io_buf),buf,var.size);
      } else if (var.nc_dtype=="int64") {
        copy_data(reinterpret_cast<long long*>(io_buf),buf,var.size);
      } else if (var.nc_dtype=="float") {
        copy_data(reinterpret_cast<float*>(io_buf),buf,var.size);
      } else if (var.nc_dtype=="double") {
        copy_data(reinterpret_cast<double*>(io_buf),buf,var.size);
      }
    }
  }
  check_scorpio_noerr (err,f.name,"variable",varname,"read_var",pioc_func);
}

// Write data from user provided buffer into the requested variable
template<typename T>
void write_var (const std::string &filename, const std::string &varname, const T* buf, const T* fillValue)
{
  EKAT_REQUIRE_MSG (buf!=nullptr,
      "Error! Cannot write in provided pointer. Invalid buffer pointer.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  const auto& f = impl::get_file(filename,"scorpio::write_var");
  auto& var = impl::get_var(filename,varname,"scorpio::write_var");

  // If the input pointer type already matches var.dtype, this is a no-op
  change_var_dtype(var,get_dtype<T>(),filename);

  int err;

  if (var.time_dep) {
    ++var.num_records;
    EKAT_REQUIRE_MSG (var.num_records==f.time_dim->length,
        "Error! Number of records for variable does not match time length.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n"
        " - time len: " + std::to_string(f.time_dim->length) + "\n"
        " - nrecords: " + std::to_string(var.num_records) + "\n");
    err = PIOc_setframe (f.ncid,var.ncid,var.num_records-1);
    check_scorpio_noerr (err,f.name,"variable",varname,"write_var","setframe");
  }

  std::string pioc_func;
  if (var.decomp) {
    // A decomposed variable, requires write_darray
    err = PIOc_write_darray(f.ncid,var.ncid,var.decomp->ncid,var.decomp->offsets.size(),buf,fillValue);
    pioc_func = "write_darray";
  } else {
    // A non-decomposed variable, use PIOc_put_var(a)
    // If nc data type doesn't match the input pointer, we need to use the var internal buffer
    const void* io_buf = buf;
    if (var.dtype!=var.nc_dtype) {
      if (var.size==-1) {
        var.size = 1;
        for (auto d : var.dims) {
          var.size *= d->length;
        }
        var.buf.resize(var.size*dtype_size(var.nc_dtype));
      }
      io_buf = var.buf.data();
      void* var_buf = var.buf.data();
      if (var.nc_dtype=="int") {
        copy_data(buf,reinterpret_cast<int*>(var_buf),var.size);
      } else if (var.nc_dtype=="int64") {
        copy_data(buf,reinterpret_cast<long long*>(var_buf),var.size);
      } else if (var.nc_dtype=="float") {
        copy_data(buf,reinterpret_cast<float*>(var_buf),var.size);
      } else if (var.nc_dtype=="double") {
        copy_data(buf,reinterpret_cast<double*>(var_buf),var.size);
      }
    }

    if (var.time_dep) {
      // We need to get the start/count for each dimension
      int ndims = var.dims.size();
      std::vector<PIO_Offset> start (ndims+1,0), count(ndims+1);
      start[0] = f.time_dim->length-1;
      count[0] = 1;
      for (int idim=0; idim<ndims; ++idim) {
        count[idim+1] = var.dims[idim]->length;
      }
      err = PIOc_put_vara(f.ncid,var.ncid,start.data(),count.data(),io_buf);
      pioc_func = "put_vara";
    } else {
      // Easy: just pass the buffer, and write all entries
      err = PIOc_put_var(f.ncid,var.ncid,io_buf);
      pioc_func = "put_var";
    }
  }
  check_scorpio_noerr (err,f.name,"variable",varname,"write_var",pioc_func);
}

// ========================== READ/WRITE ETI ========================== //

template void read_var<int>       (const std::string&, const std::string&, int*,       const int);
template void read_var<long long> (const std::string&, const std::string&, long long*, const int);
template void read_var<float>     (const std::string&, const std::string&, float*,     const int);
template void read_var<double>    (const std::string&, const std::string&, double*,    const int);
template void read_var<char>      (const std::string&, const std::string&, char*,      const int);

template void write_var<int>       (const std::string&, const std::string&, const int*,       const int*);
template void write_var<long long> (const std::string&, const std::string&, const long long*, const long long*);
template void write_var<float>     (const std::string&, const std::string&, const float*,     const float*);
template void write_var<double>    (const std::string&, const std::string&, const double*,    const double*);
template void write_var<char>      (const std::string&, const std::string&, const char*,      const char*);

// =============== Attributes operations ================== //

bool has_global_attribute (const std::string& filename, const std::string& attname)
{
  return has_attribute(filename,"GLOBAL",attname);
}

bool has_attribute (const std::string& filename, const std::string& varname, const std::string& attname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  const int ncid = pf.file->ncid;

  // Get var id
  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    const auto& var = impl::get_var(filename,varname,"scorpio::has_attribute");
    varid = var.ncid;
  }

  // Get att id
  int attid;
  int err = PIOc_inq_attid(ncid,varid,attname.c_str(),&attid);
  if (err==PIO_ENOTATT) {
    return false;
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving attribute id.\n"
      " - filename : " + filename + "\n"
      " - varname  : " + varname + "\n"
      " - attname  : " + attname + "\n"
      " - pio error: " + std::to_string(err) + "\n");

  return true;
}

template<typename T>
T get_attribute (const std::string& filename,
                 const std::string& varname,
                 const std::string& attname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = impl::get_var(filename,varname,"scorpio::get_attribute").ncid;
  }

  // If the attribute type does not match T, we need a temporary, since we can't pass T* where pio expects
  // a different type of pointer
  int att_type, err;
  err = PIOc_inq_atttype(pf.file->ncid,varid,attname.c_str(),&att_type);
  check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","inq_atttype");

  T val;
  if (att_type!=nctype(get_dtype<T>())) {

    if (att_type==PIO_INT) {
      int tmp;
      err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),reinterpret_cast<void*>(&tmp));
      check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","get_att");
      val = tmp;
    } else if (att_type==PIO_INT64) {
      std::int64_t tmp;
      err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),reinterpret_cast<void*>(&tmp));
      check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","get_att");
      val = tmp;
    } else if (att_type==PIO_FLOAT) {
      float tmp;
      err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),reinterpret_cast<void*>(&tmp));
      check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","get_att");
      val = tmp;
    } else if (att_type==PIO_DOUBLE) {
      double tmp;
      err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),reinterpret_cast<void*>(&tmp));
      check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","get_att");
      val = tmp;
    } else {
      EKAT_ERROR_MSG (
          "Unrecognized/unsupported att type\n"
          " - filename: " + filename + "\n"
          " - varname : " + varname  + "\n"
          " - attname : " + attname  + "\n"
          " - attype  : " + std::to_string(att_type) + "\n");
    }
  } else {
    err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),reinterpret_cast<void*>(&val));
    check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","get_att");
  }

  return val;
}

// Explicit instantiation
template int get_attribute (const std::string& filename,
                            const std::string& varname,
                            const std::string& attname);
template std::int64_t get_attribute (const std::string& filename,
                                     const std::string& varname,
                                     const std::string& attname);
template float get_attribute (const std::string& filename,
                              const std::string& varname,
                              const std::string& attname);
template double get_attribute (const std::string& filename,
                               const std::string& varname,
                               const std::string& attname);

// Full specialization for strings
template<>
std::string get_attribute (const std::string& filename,
                           const std::string& varname,
                           const std::string& attname)
{
  // If file wasn't open, open it on the fly. See comment in PeekFile class above.
  impl::PeekFile pf(filename);

  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = impl::get_var(filename,varname,"scorpio::set_any_attribute").ncid;
  }

  int err;
  PIO_Offset len;
  err = PIOc_inq_attlen(pf.file->ncid,varid,attname.c_str(),&len);
  check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","inq_attlen");

  std::string val(len,'\0');

  err = PIOc_get_att(pf.file->ncid,varid,attname.c_str(),val.data());
  check_scorpio_noerr(err,filename,"attribute",attname,"get_attribute","put_att");

  return val;
}

template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const T& att)
{
  const auto& f = impl::get_file (filename,"scorpio::set_any_attribute");

  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    varid = impl::get_var(filename,varname,"scorpio::set_any_attribute").ncid;
  }

  // If the file was not in define mode, we must call enddef at the end
  const bool needs_redef = f.enddef;
  if (needs_redef) {
    redef(filename);
  }
  
  int err = PIOc_put_att(f.ncid,varid,attname.c_str(),nctype<T>(),nclen(att),ncdata(att));
  check_scorpio_noerr(err,filename,"attribute",attname,"set_attribute","put_att");

  if (needs_redef) {
    enddef(filename);
  }
}

// Explicit instantiation
template void set_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname,
                             const int& att);
template void set_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname,
                             const std::int64_t& att);
template void set_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname,
                             const float& att);
template void set_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname,
                             const double& att);
template void set_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname,
                             const std::string& att);

} // namespace scorpio
} // namespace scream
