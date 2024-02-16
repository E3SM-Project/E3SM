#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "share/util/scream_time_stamp.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <string>

namespace scream
{

enum class FileType {
  ModelOutput,
  ModelRestart,
  HistoryRestart,
  Unset
};

inline std::string e2str(const FileType avg) {
  using FT = FileType;
  switch (avg) {
    case FT::ModelOutput:     return "model-output";
    case FT::ModelRestart:    return "model-restart";
    case FT::HistoryRestart:  return "history-restart";
    default:                  return "UNSET";
  }
}

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

  // If frequency_units is not "none" or "never", frequency *must* be set to a positive number
  int frequency = -1;
  std::string frequency_units = "none";

  int nsamples_since_last_write = 0;  // Needed when updating output data, such as with the OAT::Average flag

  util::TimeStamp next_write_ts;
  util::TimeStamp last_write_ts;

  mutable int dt = 0;

  bool output_enabled () const {
    return frequency_units!="none" && frequency_units!="never";
  }

  void compute_dt (const util::TimeStamp& ts) {
    if (dt==0 and frequency_units=="nsteps") {
      int nsteps = ts.get_num_steps() - last_write_ts.get_num_steps();
      if (nsteps>0) {
        dt = (ts-last_write_ts) / nsteps;
      }

      // Now that we have dt, we can compute the next write ts
      // NOTE: when OutputManager calls this for t0 output, we still have dt=0,
      //       so this call sets next_write_ts=last_write_ts. We'll have to wait
      //       until the first call during the time loop for a nonzero dt.
      compute_next_write_ts();
    }
  }

  bool is_write_step (const util::TimeStamp& ts) const {
    if (not output_enabled()) return false;
    return ts==next_write_ts;
  }

  // Computes next_write_ts from frequency and last_write_ts
  void compute_next_write_ts () {
    EKAT_REQUIRE_MSG (last_write_ts.is_valid(),
        "Error! Cannot compute next_write_ts, since last_write_ts was never set.\n");
    if (frequency_units=="nsteps") {
      next_write_ts = last_write_ts;
      if (dt>0) {
        next_write_ts += dt*frequency;
        next_write_ts.set_num_steps(last_write_ts.get_num_steps()+frequency);
      }
    } else if (frequency_units=="nsecs") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency;
    } else if (frequency_units=="nmins") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*60;
    } else if (frequency_units=="nhours") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*3600;
    } else if (frequency_units=="ndays") {
      next_write_ts = last_write_ts;
      next_write_ts += frequency*86400;
    } else if (frequency_units=="nmonths" or frequency_units=="nyears") {
      auto date = last_write_ts.get_date();
      if (frequency_units=="nmonths") {
        int temp = date[1] + frequency - 1;
        date[1] = temp % 12 + 1;
        date[0] += temp / 12;
      } else {
        date[0] += frequency;
      }

      // Fix day, in case we moved to a month/year where current days. E.g., if last_write
      // was on Mar 31st, and units='nmonths', next write is on Apr 30th. HOWEVER, this
      // means we will *always* write on the 30th of each month after then, since we have
      // no memory of the fact that we were writing on the 31st before.
      date[2] = std::min(date[2],util::days_in_month(date[0],date[1]));
      next_write_ts = util::TimeStamp(date,last_write_ts.get_time());
    } else {
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported frequency unit '" + frequency_units + "'\n");
    }
  }

  void update_write_timestamps () {
    // Shift next_write_ts into last_write_ts, recompute next_write_ts, and zero-out nsamples_since_last_write
    last_write_ts = next_write_ts;
    compute_next_write_ts ();
    nsamples_since_last_write = 0;
  }
};

// Mini struct to hold some specs of an IO file
// To keep nc files small, we limit the number of snapshots in each nc file
// When the number of snapshots in a file reaches m_out_max_steps, it's time
// to close the out file, and open a new one.
struct IOFileSpecs {
  bool is_open = false;
  std::string filename;
  int num_snapshots_in_file = 0;
  int max_snapshots_in_file;

  // If positive, flush the output file every these many snapshots
  int flush_frequency = -1;

  bool file_is_full () const { return num_snapshots_in_file>=max_snapshots_in_file; }
  bool file_needs_flush () const { return flush_frequency>0 and num_snapshots_in_file%flush_frequency==0; }

  // Adding number of MPI ranks to the filenamea is useful in testing, since we can run
  // multiple instances of the same test in parallel (with different number of ranks),
  // without the risk of them overwriting each other output.
  // For production runs, this is not desirable.
  bool filename_with_mpiranks    = false;

  bool save_grid_data            = true;

  // Whether it is a model output, model restart, or history restart file
  FileType ftype = FileType::Unset;

  bool is_restart_file () const {
    return ftype==FileType::ModelRestart or ftype==FileType::HistoryRestart;
  }

  std::string suffix () const {
    if (ftype==FileType::HistoryRestart)
      return ".rhist";
    else if (ftype==FileType::ModelRestart)
      return ".r";
    else
      return "";
  }
};

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
