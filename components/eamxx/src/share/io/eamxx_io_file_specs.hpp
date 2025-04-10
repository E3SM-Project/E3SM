#ifndef SCREAM_IO_FILE_SPECS_HPP
#define SCREAM_IO_FILE_SPECS_HPP

#include "share/io/eamxx_scorpio_types.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>
#include <limits>

namespace scream
{

// How the file capacity is specified
// NOTE: Keep Yearly=0, Monthly=1, Daily=2, so you can use them
//       to access a TimeStamp date at the correct index in StorageSpecs methods
enum StorageType {
  Yearly  = 0,  // Each file contains output for one year
  Monthly = 1,  // Each file contains output for one month
  Daily   = 2,  // Each file contains output for one day
  NumSnaps      // Fixed number of snaps per file
};

inline std::string e2str (const StorageType st) {
  switch (st) {
    case NumSnaps: return "num_snapshots";
    case Yearly:   return "one_year";
    case Monthly:  return "one_month";
    case Daily:    return "one_day";
    default:       return "unknown";
  }
}

struct StorageSpecs {

  StorageType type = NumSnaps;

  // Current index for type!=NumSnaps. It stores the year/month/day
  // index associated with this file
  int  curr_idx = -1;

  // A snapshot fits if
  //  - type=NumSnaps: the number of stored snaps is less than the max allowed per file.
  //  - otherwise: the snapshot year/month/day index match the one currently stored in the file
  //               or the file has no snapshot stored yet
  bool snapshot_fits (const util::TimeStamp& t) const {
    switch (type) {
      case Yearly:  [[fallthrough]];
      case Monthly: [[fallthrough]];
      case Daily:
        return curr_idx==-1 or curr_idx==t.get_date()[static_cast<int>(type)];
      case NumSnaps:
        return num_snapshots_in_file<max_snapshots_in_file;
      default:
        EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
    }
    return false;
  }

  void update_storage (const util::TimeStamp& t) {
    // We always update the snap counter (regardless of storage type),
    // so that FileSpecs can correctly detect if it needs flushing
    switch (type) {
      case Yearly:  [[fallthrough]];
      case Monthly: [[fallthrough]];
      case Daily:
          curr_idx = t.get_date()[static_cast<int>(type)]; [[fallthrough]];
      case NumSnaps:
        ++num_snapshots_in_file;  break;
      default:
        EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
    }
  }

  // If type==NumSnaps, we need to keep track of the snapshots count
  int num_snapshots_in_file =  0;
  int max_snapshots_in_file = std::numeric_limits<int>::max();
};

// Mini struct to hold some specs of an IO file
// To keep nc files small, we limit the number of snapshots in each nc file
// When the number of snapshots in a file reaches m_out_max_steps, it's time
// to close the out file, and open a new one.
struct IOFileSpecs {

  StorageSpecs storage = {};

  bool is_open = false;
  std::string filename;

  scorpio::IOType iotype = scorpio::IOType::Invalid;

  // If positive, flush the output file every these many snapshots
  int flush_frequency = std::numeric_limits<int>::max();

  bool file_needs_flush () const {
    return storage.num_snapshots_in_file>0 and storage.num_snapshots_in_file%flush_frequency==0;
  }

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

  void close () {
    is_open = false;
    storage.num_snapshots_in_file = 0;
    storage.curr_idx = -1;
  }
};

} // namespace scream
#endif // SCREAM_IO_FILE_SPECS_HPP
