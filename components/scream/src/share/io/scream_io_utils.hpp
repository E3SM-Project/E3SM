#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "share/util/scream_time_stamp.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include <string>

namespace scream
{

enum class OutputAvgType {
  Instant,
  Max,
  Min,
  Average,
  Invalid
};

inline std::string e2str(const OutputAvgType avg) {
  using OAT = OutputAvgType;
  switch (avg) {
    case OAT::Instant:  return "INSTANT";
    case OAT::Max:      return "MAX";
    case OAT::Min:      return "MIN";
    case OAT::Average:  return "AVERAGE";
    default:            return "INVALID";
  }
}

inline OutputAvgType str2avg (const std::string& s) {
  auto s_ci = ekat::upper_case(s);
  using OAT = OutputAvgType;
  for (auto e : {OAT::Instant, OAT::Max, OAT::Min, OAT::Average}) {
    if (s_ci==e2str(e)) {
      return e;
    }
  }

  return OAT::Invalid;
}

// Mini struct to hold IO frequency info
struct IOControl {
  // A non-positive frequency can be used to signal IO disabled
  int frequency = -1;
  int nsteps_since_last_write;
  std::string frequency_units;

  bool is_write_step () {
    return frequency>0 && (nsteps_since_last_write % frequency == 0);
  }
};

// Mini struct to hold some specs of an IO file
// To keep nc files small, we limit the number of snapshots in each nc file
// When the number of snapshots in a file reaches m_out_max_steps, it's time
// to close the out file, and open a new one.
struct IOFileSpecs {
  bool is_open = false;
  std::string filename;
  int num_snapshots_in_file;
  int max_snapshots_in_file;
  bool file_is_full () const { return num_snapshots_in_file==max_snapshots_in_file; }
  // Whether a time string or the number of mpiranks should be attached to the filename.
  // Attaching the timestamp to each file ensures that in runs with multiple
  // outputs, old output isn't accidentally overwritten, so it's on by default.
  // Attaching the nranks, OTOH, is necessary only when running multiple mpi configs
  // of the same test, so they run concurrently, writing on different files.
  bool filename_with_time_string = true;
  bool filename_with_mpiranks    = false;
  bool filename_with_avg_type    = true;
  bool filename_with_frequency   = true;
};

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
