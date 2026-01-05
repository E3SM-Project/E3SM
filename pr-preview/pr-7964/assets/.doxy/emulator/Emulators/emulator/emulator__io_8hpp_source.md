

# File emulator\_io.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_io.hpp**](emulator__io_8hpp.md)

[Go to the documentation of this file](emulator__io_8hpp.md)


```C++


#ifndef EMULATOR_IO_HPP
#define EMULATOR_IO_HPP

#include <mpi.h>
#include <string>
#include <vector>

// PIO C interface from SCORPIO
extern "C" {
#include <pio.h>
}

namespace emulator {

class EmulatorIO {
public:
  // =========================================================================
  // Initialization
  // =========================================================================

  static void initialize(MPI_Comm comm, const std::string &comp_name);

  static void finalize();

  static bool is_initialized() { return s_initialized; }

  // =========================================================================
  // File Operations
  // =========================================================================

  static int open_file(const std::string &filename);

  static int create_file(const std::string &filename);

  static void close_file(int ncid);

  static void sync_file(int ncid);

  // =========================================================================
  // Variable I/O
  // =========================================================================

  static bool read_var_1d(int ncid, const std::string &varname, double *data,
                          int size);

  static bool read_var_2d(int ncid, const std::string &varname, double *data,
                          int nx, int ny);

  static bool read_var_3d_slice(int ncid, const std::string &varname,
                                double *data, int nx, int ny, int time_idx);

  static bool read_var_1d_int(int ncid, const std::string &varname, int *data,
                              int size);

  static bool write_var_1d(int ncid, const std::string &varname,
                           const double *data, int size);

  static bool write_var_2d(int ncid, const std::string &varname,
                           const double *data, int nx, int ny);

  // =========================================================================
  // Dimension/Variable Operations
  // =========================================================================

  static int define_dim(int ncid, const std::string &dimname, int length);

  static int get_dim_size(int ncid, const std::string &dimname);

  static bool has_var(int ncid, const std::string &varname);

  static int define_var(int ncid, const std::string &varname, int nctype,
                        const std::vector<int> &dimids);

  static bool end_def(int ncid);

private:
  static MPI_Comm s_comm;    
  static int s_iosysid;      
  static bool s_initialized; 
  static int s_rank;         
};

constexpr double FILLVALUE = 1.0e20;

} // namespace emulator

#endif // EMULATOR_IO_HPP
```


