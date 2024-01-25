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

inline int str2iotype(const std::string &str)
{
  if(str == "default"){
    return 0;
  }
  else if(str == "netcdf"){
    return 1;
  }
  else if(str == "pnetcdf"){
    return 2;
  }
  else if(str == "adios"){
    return 3;
  }
  else if(str == "hdf5"){
    return 4;
  }
  else{
    return 0;
  }
}

inline std::string iotype2str(int iotype)
{
  switch(iotype){
    case 0: return "default";
    case 1: return "netcdf";
    case 2: return "pnetcdf";
    case 3: return "adios";
    case 4: return "hdf5";
    default: return "default";
  }
}

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
