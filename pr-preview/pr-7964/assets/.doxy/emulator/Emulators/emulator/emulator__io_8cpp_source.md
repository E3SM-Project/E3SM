

# File emulator\_io.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_io.cpp**](emulator__io_8cpp.md)

[Go to the documentation of this file](emulator__io_8cpp.md)


```C++


#include "emulator_io.hpp"

namespace emulator {

// Static member definitions
MPI_Comm EmulatorIO::s_comm = MPI_COMM_NULL;
int EmulatorIO::s_iosysid = -1;
bool EmulatorIO::s_initialized = false;
int EmulatorIO::s_rank = -1;

void EmulatorIO::initialize(MPI_Comm comm, const std::string &comp_name) {
  (void)comp_name; // Currently unused

  if (s_initialized) {
    return;
  }

  s_comm = comm;

  // Check if MPI communicator is valid
  int nprocs = 1;
  if (comm != MPI_COMM_NULL && comm != 0) {
    MPI_Comm_rank(comm, &s_rank);
    MPI_Comm_size(comm, &nprocs);
  } else {
    // Fallback for non-MPI runs or invalid comm
    s_rank = 0;
    return; // Cannot initialize PIO without valid MPI
  }

  // PIO initialization parameters
  int stride = 1;
  int base = 0;
  int rearranger = PIO_REARR_SUBSET;

  int ret =
      PIOc_Init_Intracomm(comm, nprocs, stride, base, rearranger, &s_iosysid);
  if (ret != PIO_NOERR) {
    // Initialization failed - caller should handle gracefully
    return;
  }

  s_initialized = true;
}

void EmulatorIO::finalize() {
  if (!s_initialized) {
    return;
  }

  if (s_iosysid >= 0) {
    PIOc_finalize(s_iosysid);
    s_iosysid = -1;
  }

  s_initialized = false;
}

int EmulatorIO::open_file(const std::string &filename) {
  if (!s_initialized) {
    return -1;
  }

  int ncid;
  int iotype = PIO_IOTYPE_NETCDF;
  int ret =
      PIOc_openfile(s_iosysid, &ncid, &iotype, filename.c_str(), PIO_NOWRITE);
  if (ret != PIO_NOERR) {
    return -1;
  }
  return ncid;
}

int EmulatorIO::create_file(const std::string &filename) {
  if (!s_initialized) {
    return -1;
  }

  int ncid;
  int iotype = PIO_IOTYPE_NETCDF;
  int ret =
      PIOc_createfile(s_iosysid, &ncid, &iotype, filename.c_str(), PIO_CLOBBER);
  if (ret != PIO_NOERR) {
    return -1;
  }
  return ncid;
}

void EmulatorIO::close_file(int ncid) {
  if (ncid >= 0) {
    PIOc_closefile(ncid);
  }
}

void EmulatorIO::sync_file(int ncid) {
  if (ncid >= 0) {
    PIOc_sync(ncid);
  }
}

bool EmulatorIO::read_var_1d(int ncid, const std::string &varname, double *data,
                             int size) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  PIO_Offset start[1] = {0};
  PIO_Offset count[1] = {static_cast<PIO_Offset>(size)};

  ret = PIOc_get_vara_double(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

bool EmulatorIO::read_var_2d(int ncid, const std::string &varname, double *data,
                             int nx, int ny) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  PIO_Offset start[2] = {0, 0};
  PIO_Offset count[2] = {static_cast<PIO_Offset>(ny),
                         static_cast<PIO_Offset>(nx)};

  ret = PIOc_get_vara_double(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

bool EmulatorIO::read_var_3d_slice(int ncid, const std::string &varname,
                                   double *data, int nx, int ny, int time_idx) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  // Read slice at time_idx from (time, lat, lon) array
  PIO_Offset start[3] = {static_cast<PIO_Offset>(time_idx), 0, 0};
  PIO_Offset count[3] = {1, static_cast<PIO_Offset>(ny),
                         static_cast<PIO_Offset>(nx)};

  ret = PIOc_get_vara_double(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

bool EmulatorIO::read_var_1d_int(int ncid, const std::string &varname,
                                 int *data, int size) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  PIO_Offset start[1] = {0};
  PIO_Offset count[1] = {static_cast<PIO_Offset>(size)};

  ret = PIOc_get_vara_int(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

bool EmulatorIO::write_var_1d(int ncid, const std::string &varname,
                              const double *data, int size) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  PIO_Offset start[1] = {0};
  PIO_Offset count[1] = {static_cast<PIO_Offset>(size)};

  ret = PIOc_put_vara_double(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

bool EmulatorIO::write_var_2d(int ncid, const std::string &varname,
                              const double *data, int nx, int ny) {
  int varid;
  int ret = PIOc_inq_varid(ncid, varname.c_str(), &varid);
  if (ret != PIO_NOERR) {
    return false;
  }

  PIO_Offset start[2] = {0, 0};
  PIO_Offset count[2] = {static_cast<PIO_Offset>(ny),
                         static_cast<PIO_Offset>(nx)};

  ret = PIOc_put_vara_double(ncid, varid, start, count, data);
  return (ret == PIO_NOERR);
}

int EmulatorIO::define_dim(int ncid, const std::string &dimname, int length) {
  int dimid;
  int ret = PIOc_def_dim(ncid, dimname.c_str(), length, &dimid);
  return (ret == PIO_NOERR) ? dimid : -1;
}

int EmulatorIO::get_dim_size(int ncid, const std::string &dimname) {
  int dimid;
  int ret = PIOc_inq_dimid(ncid, dimname.c_str(), &dimid);
  if (ret != PIO_NOERR) {
    return -1;
  }

  PIO_Offset dimlen;
  ret = PIOc_inq_dimlen(ncid, dimid, &dimlen);
  return (ret == PIO_NOERR) ? static_cast<int>(dimlen) : -1;
}

bool EmulatorIO::has_var(int ncid, const std::string &varname) {
  int varid;
  return (PIOc_inq_varid(ncid, varname.c_str(), &varid) == PIO_NOERR);
}

int EmulatorIO::define_var(int ncid, const std::string &varname, int nctype,
                           const std::vector<int> &dimids) {
  int varid;
  int ret =
      PIOc_def_var(ncid, varname.c_str(), nctype,
                   static_cast<int>(dimids.size()), dimids.data(), &varid);
  return (ret == PIO_NOERR) ? varid : -1;
}

bool EmulatorIO::end_def(int ncid) {
  int ret = PIOc_enddef(ncid);
  return (ret == PIO_NOERR);
}

} // namespace emulator
```


