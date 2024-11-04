#ifndef SCREAM_IO_FILE_SPECS_HPP
#define SCREAM_IO_FILE_SPECS_HPP

#include "share/io/scream_scorpio_types.hpp"
#include "share/io/scream_io_utils.hpp"
#include "share/util/scream_time_stamp.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>
#include <limits>

namespace scream
{

// How the file capacity is specified
enum StorageType {
  NumSnaps,   // Fixed number of snaps per file
  Monthly,    // Each file contains output for one month
  Yearly      // Each file contains output for one year
};

inline std::string e2str (const StorageType st) {
  switch (st) {
    case NumSnaps: return "num_snapshots";
    case Yearly:   return "one_year";
    case Monthly:  return "one_month";
    default:       return "unknown";
  }
}

struct StorageSpecs {

  StorageType type = NumSnaps;

  // Current index ***in terms of this->type***
  // If type==NumSnaps, curr_idx=num_snapshots_in_file,
  // otherwise it is the month/year index stored in this file
  int  curr_idx = -1;

  // A snapshot fits if
  //  - type=NumSnaps: the number of stored snaps is less than the max allowed per file.
  //  - otherwise: the snapshot month/year index match the one currently stored in the file
  //               or the file has no snapshot stored yet
  bool snapshot_fits (const util::TimeStamp& t) {
    const auto& idx = type==Monthly ? t.get_month() : t.get_year();
    switch (type) {
      case Yearly:
      case Monthly:
        return curr_idx==-1 or curr_idx==idx;
      case NumSnaps:
        return num_snapshots_in_file<max_snapshots_in_file;
      default:
        EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
    }
  }

  void update_storage (const util::TimeStamp& t) {
    switch (type) {
      case Yearly:    curr_idx = t.get_year();  break;
      case Monthly:   curr_idx = t.get_month(); break;
      case NumSnaps:  ++num_snapshots_in_file;  break;
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

  // bool file_is_full () const { return num_snapshots_in_file>=max_snapshots_in_file; }
  bool file_needs_flush () const {
    return storage.num_snapshots_in_file%flush_frequency==0;
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
