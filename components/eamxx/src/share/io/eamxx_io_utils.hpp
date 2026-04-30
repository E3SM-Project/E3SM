#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "eamxx_io_control.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/grid/abstract_grid.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_comm.hpp>

#include <string>
#include <memory>
#include <utility>

namespace scream
{

// Universal fill value for variables
// NOTE: for floating point numbers, use the SAME numerical value, so that
//       we don't need to be aware of the precision of variables when checking
//       against fill_value (e.g., when reading in double data that was saved
//       in single precision)
template<typename T>
constexpr T fill_value =
  std::is_integral_v<T> ? std::numeric_limits<int>::max() / 2
                        : std::is_floating_point_v<T> ? std::numeric_limits<float>::max() / static_cast<float>(1e5)
                                                      : std::numeric_limits<char>::max();

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

// The AD will pass a default constructed control, since it doesn't know the values
// of REST_N/REST_OPTION used in the previous run
// Output streams MUST pass a valid control structure, cause we need to differentiate
// between, e.g., streams with same filename prefix, but different output freq specs
std::string find_filename_in_rpointer (
    const std::string& filename_prefix,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0,
    const bool allow_not_found = false,
    const OutputAvgType avg_type = OutputAvgType::Instant,
    const IOControl& control = {});

// Shortcut to write/read to/from YYYYMMDD/HHMMSS attributes in the NC file
void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps = false);
util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps = false);

// Parse a CF-compliant time units string of the form "<unit> since <date> [<time>]"
// and return the reference TimeStamp and multiplier (in seconds) for the given unit.
// Supported units: seconds, minutes, hours, days.
// Supported date formats: YYYY-MM-DD, YYYY-MM-DD HH:MM:SS, YYYY-MM-DD-SSSSS
// NOTE: the filename is only used to give context if throwing an exception
std::pair<util::TimeStamp,int>
parse_cf_time_units (const std::string& units_str,
                     const std::string& filename);

// Create a diagnostic from a string representation of it.
// E.g., create the diag to compute fieldX_at_500hPa.
std::shared_ptr<AtmosphereDiagnostic>
create_diagnostic (const std::string& diag_name,
                   const std::shared_ptr<const AbstractGrid>& grid);

} // namespace scream

#endif // SCREAM_IO_UTILS_HPP
