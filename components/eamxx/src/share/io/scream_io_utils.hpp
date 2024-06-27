#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "share/util/scream_time_stamp.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/mpi/ekat_comm.hpp>

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

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0);

struct LongNames {

  std::string get_longname (const std::string& name) {
    if (name_2_longname.count(name)>0) {
      return name_2_longname.at(name);
    } else {
      // TODO: Do we want to print a Warning message?  I'm not sure if its needed.
      return name;
    }
  }

  // Create map of longnames, can be added to as developers see fit.
  std::map<std::string,std::string> name_2_longname = {
	  {"lev","hybrid level at midpoints (1000*(A+B))"},
	  {"hyai","hybrid A coefficient at layer interfaces"},
    {"hybi","hybrid B coefficient at layer interfaces"},
    {"hyam","hybrid A coefficient at layer midpoints"},
    {"hybm","hybrid B coefficient at layer midpoints"}
  };
  
};

// Shortcut to write/read to/from YYYYMMDD/HHMMSS attributes in the NC file
void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps = false);
util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps = false);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
