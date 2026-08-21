#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "eamxx_io_control.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/grid/abstract_grid.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_comm.hpp>

#include <string>
#include <memory>
#include <utility>
#include <vector>

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

// Create a diagnostic from a string representation of it.
// E.g., create the diag to compute fieldX_at_500hPa.
std::shared_ptr<AbstractDiagnostic>
create_diagnostic (const std::string& diag_name,
                   const std::shared_ptr<const AbstractGrid>& grid);

// Copy every parameter of 'src' into 'dst', overwriting what is already there.
// Used to overlay the per-diag options given in an output yaml on top of the
// params a diagnostic expression lowered to. Only the scalar and vector types
// the yaml parser can produce are handled; anything else throws, rather than
// being silently dropped.
void overlay_params (ekat::ParameterList& dst, const ekat::ParameterList& src);

// Same, but only for the parameters named in 'only'. A name that 'src' does not
// have is skipped.
void overlay_params (ekat::ParameterList& dst, const ekat::ParameterList& src,
                     const std::vector<std::string>& only);

// A deterministic string rendering of 'params', for comparing two parameter
// lists for equality. Keys are sorted, so the result does not depend on the
// order the yaml listed them in.
std::string params_signature (const ekat::ParameterList& params);

} // namespace scream

#endif // SCREAM_IO_UTILS_HPP
