#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "eamxx_scorpio_types.hpp"

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/ekat_assert.hpp>

#include <string>
#include <vector>

/*
 * This file contains interfaces to scorpio C library routines
 * 
 * There are two reasons for these interfaces:
 * 
 * - allow better context when exceptions/errors occour: throwing exceptions
 *   can even allow host code to use try-catch blocks to figure out what works.
 *   Moreover, the user does not need to continuously check return codes, since
 *   we throw any time scorpio returns something different from PIO_NOERR.
 * - allow using simpler interfaces: by using string/vector and templates,
 *   the host code can have more compact interfaces.
 *
 * Most of the "get" interfaces can open a file in read mode on the fly.
 * This allows easier usage from the user point of view, but has some drawbacks:
 *  - the file is open/closed every time you query it
 *  - the file is open with default IOType
 * If you need to do several queries, you should consider registering the file manually
 * before doing any query, and release it once you are done.
 */

namespace scream {
namespace scorpio {

inline std::string default_time_name () { return "time"; }

// Handles aliases for same type:
//  single -> float
//  float  -> float
//  double -> double
//  real   -> float or double (depending on SCREAM_DOUBLE_PRECISION)
std::string refine_dtype (const std::string& dtype);

template<typename T>
std::string get_dtype () {
  using raw_t = typename std::remove_cv<T>::type;
  std::string s;
  if constexpr (std::is_same<raw_t,int>::value) {
    s = "int";
  } else if constexpr (std::is_same<raw_t,float>::value) {
    s = "float";
  } else if constexpr (std::is_same<raw_t,double>::value) {
    s = "double";
  } else if constexpr (std::is_integral<raw_t>::value &&
              std::is_signed<raw_t>::value &&
              sizeof(raw_t)==sizeof(long long)) {
    s = "int64";
  } else if constexpr (std::is_same<raw_t,char>::value) {
    s = "char";
  } else {
    EKAT_ERROR_MSG ("Error! Invalid/unsupported data type.\n");
  }
  return s;
}

// =================== Global operations ================= //

void init_subsystem(const ekat::Comm& comm, const int atm_id = 0);
bool is_subsystem_inited ();
void finalize_subsystem ();

// =================== File operations ================= //

// Opens a file, returns const handle to it (useful for Read mode, to get dims/vars)
void register_file(const std::string& filename, const FileMode mode, const IOType iotype = IOType::DefaultIOType);

// Release a file (if in Write mode, sync and close the file);
void release_file  (const std::string& filename);

// Check if file is open. If mode!=Unset, also checks that it's open with given mode
bool is_file_open (const std::string& filename, const FileMode mode = Unset);

// Force a flush to file (for Write mode only)
void flush_file (const std::string &filename);

// Reopen/ends the definition phase
void redef (const std::string &filename);
void enddef (const std::string &filename);

// =================== Dimensions operations ======================= //

// Define dim on output file (cannot call on Read/Append files)
void define_dim (const std::string& filename, const std::string& dimname, const int length);

// Check that the given dimension is in the file. If length>0, also check that the length is as expected.
bool has_dim (const std::string& filename,
              const std::string& dimname,
              const int length = -1);

int get_dimlen (const std::string& filename, const std::string& dimname);
int get_dimlen_local (const std::string& filename, const std::string& dimname);

// When we read/write a var, we may want to only process one slice of the var,
// orresponding to a particular time index. In this case, the time dim is a
// "record" dimension. When we open a file, if there is an unlimited dim, it is
// automatically marked as the "time" dim (regardless of its name). But if there
// is no unlimited dim, we may need to tell the interface to interpret one of the
// dims as the record dimension (which in the interface we refer to as "time" dim).
void mark_dim_as_time (const std::string& filename, const std::string& dimname);
bool has_time_dim (const std::string& filename);

// This is used by I/O when restarting a simulation after a crash: the existing file
// may contain some snapshots after the rest time, but we may want to overwrite them.
void reset_time_dim_len(const std::string& filename, const int new_length);

// Get len/name of the time dimension
// NOTE: these throw if time dim is not present. Use has_dim to check first.
int get_time_len (const std::string& filename);
std::string get_time_name (const std::string& filename);

// =================== Decompositions operations ==================== //

// Create a decomposition along a particular dimension
// Notes:
// - we declare a decomposition along a single dimension.
// - set_dim_decomp requires *offsets* in the global array, not global indices
//   (for 0-based indices, they're the same, but for 1-based indices they're not)
// - the second version is a shortcut for contiguous decompositions. Notice that NO
//   check is performed to ensure that the decomposition covers the entire dimension,
//   or that there are no overlaps (in fact, one *may* need overlap, in certain reads).
// - the third version is a shortcut of the second, where we compute start/count based
//   on a linear decomposition of the dimension along all ranks in the IO comm stored
//   in the ScorpioInstance. The return value is the local length of the dimension
// - if allow_reset=true, we simply reset the decomposition (if present).
// - if allow_reset=false, if a decomposition for this dim is already set, we error out

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const std::vector<offset_t>& my_offsets,
                     const bool allow_reset = false);

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const offset_t start, const offset_t count,
                     const bool allow_reset = false);

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const bool allow_reset = false);

// ================== Variable operations ================== //

// Define var on output file (cannot call on Read/Append files)
void define_var (const std::string& filename, const std::string& varname,
                 const std::string& units, const std::vector<std::string>& dimensions,
                 const std::string& dtype, const std::string& nc_dtype,
                 const bool time_dependent = false);

// Shortcut when units are not used, and dtype==nc_dtype
void define_var (const std::string& filename, const std::string& varname,
                 const std::vector<std::string>& dimensions,
                 const std::string& dtype,
                 const bool time_dependent = false);

// This is useful when reading data sets. E.g., if the pio file is storing
// a var as float, but we need to read it as double, we need to call this.
// NOTE: read_var/write_var automatically change the dtype if the input
//       pointer type does not match the var dtype. However, changing dtype
//       forces a rebuild of the var decomp (if any). Hence, if you know
//       the var WILL be read/written as decomposed, you should call this method
//       BEFORE calling set_dim_decomp, so that the decomp is built directly
//       with the correct data type (PIO decomps depend on var dtype).
void change_var_dtype (const std::string& filename,
                       const std::string& varname,
                       const std::string& dtype);

// Check that the given variable is in the file.
bool has_var (const std::string& filename, const std::string& varname);

// Allows to easily query var metadata, such as dims, units, data type,...
const PIOVar& get_var (const std::string& filename,
                       const std::string& varname);

// Defines both a time dimension and a time variable
void define_time (const std::string& filename, const std::string& units,
                  const std::string& time_name = default_time_name());

// Update value of time variable, increasing time dim length
void update_time(const std::string &filename, const double time);

// Retrieves the time variable value(s)
double get_time (const std::string& filename, const int time_index = -1);
std::vector<double> get_all_times (const std::string& filename);

// Read variable into user provided buffer.
// If time dim is present, read given time slice (time_index=-1 means "read last record).
// If time dim is not present and time_index>=0, it is interpreted as the index of the
// first dimension (which is not unlimited).
// NOTE: ETI in the cpp file for int, float, double.
template<typename T>
void read_var (const std::string &filename, const std::string &varname, T* buf, const int time_index = -1);

// Write data from user provided buffer into the requested variable
// NOTE: ETI in the cpp file for int, float, double.
template<typename T>
void write_var (const std::string &filename, const std::string &varname, const T* buf, const T* fillValue = nullptr);

// =============== Attributes operations ================== //

// To specify GLOBAL attributes, pass "GLOBAL" as varname
// NOTE: set/get_attribute are implemented in the cpp file,
//       with Explicit Instantiation only for:
//         int, std::int64_t, float, double, std::string

bool has_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname);

template<typename T>
T get_attribute (const std::string& filename,
                 const std::string& varname,
                 const std::string& attname);

template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const T& att);

// Shortcut, to allow calling set_attribute with compile-time strings, like so
//   set_attribute(my_file,my_var,my_att_name,"my_value");
template<int N>
inline void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const char (&att)[N])
{
  set_attribute<std::string>(filename,varname,attname,att);
}

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
