#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "share/util/scream_time_stamp.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include <string>
#include "share/util/scream_time_stamp.hpp"

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
  int nsamples_since_last_write;  // Needed when updating output data, such as with the OAT::Average flag
  util::TimeStamp timestamp_of_last_write;
  std::string frequency_units = "none";

  bool is_write_step (const util::TimeStamp& ts) {
    // Mini-routine to determine if it is time to write output to file.
    // The current allowable options are nsteps, nsecs, nmins, nhours, ndays, nmonths, nyears
    // We query the frequency_units string value to determine which option it is.
    bool ret = false;
    if (frequency > 0 && frequency_units != "never" && frequency_units != "none") {
      auto ts_diff = (ts-timestamp_of_last_write);
      if (frequency_units == "nsteps") {
        // Just use the num_steps from timestamps
        return ((ts.get_num_steps()-timestamp_of_last_write.get_num_steps()) % frequency == 0);
      // We will need to use timestamp information
      } else if (frequency_units == "nsecs") {
        ret = ((ts_diff > 0) && (ts_diff % frequency == 0));
      } else if (frequency_units == "nmins") {
        ret = (ts_diff >= 60) && (ts_diff % (frequency*60) == 0);
      } else if (frequency_units == "nhours") {
        ret = (ts_diff >= 3600) && (ts_diff % (frequency*3600) == 0);
      } else if (frequency_units == "ndays") {
        ret = (ts_diff >= 86400) && (ts_diff % (frequency*86400) == 0);
      } else if (frequency_units == "nmonths" || frequency_units == "nyears") {
        // For months and years we need to be careful, can't just divide ts_diff by a set value.
        // First we make sure that if we are the same day of the month and at the same time of day.
        if (ts.get_day() == timestamp_of_last_write.get_day() &&
            ts.sec_of_day() == timestamp_of_last_write.sec_of_day()) {
          auto diff = 0;
          // Determine how many years have passed
          diff += (ts.get_year()  - timestamp_of_last_write.get_year());
          if  (frequency_units == "nyears") {
            ret = (diff>0) && (diff % frequency == 0);
          }
          // Determine number of months that have passed
          diff *= 12;
          diff += (ts.get_month() - timestamp_of_last_write.get_month());
          if  (frequency_units == "nmonths") {
            ret = (diff>0) && (diff % frequency == 0);
          }
        } 
      } else {
        EKAT_REQUIRE_MSG(false,"Invalid frequency unit of [" + frequency_units + "] for output stream.  Please check that all outputs have frequency_units of\n"
                               "none, never, nsteps, nsecs, nmins, nhours, ndays, nmonths, nyears");
      }
    }
    return ret;
  } // End function is_write_step
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
  bool save_grid_data            = true;
};

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
