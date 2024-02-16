#ifndef SCREAM_IO_FILE_SPECS_HPP
#define SCREAM_IO_FILE_SPECS_HPP

#include "share/io/scream_io_utils.hpp"
#include "share/util/scream_time_stamp.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>
#include <limits>

namespace scream
{

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

} // namespace scream
#endif // SCREAM_IO_FILE_SPECS_HPP
