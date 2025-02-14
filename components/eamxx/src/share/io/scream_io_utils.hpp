#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "scream_io_control.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/grid/abstract_grid.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <string>
#include <memory>
#include <regex>

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
    const OutputAvgType avg_type = OutputAvgType::Instant,
    const IOControl& control = {});

// Shortcut to write/read to/from YYYYMMDD/HHMMSS attributes in the NC file
void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps = false);
util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps = false);

// A method to set diagnostics (mainly to organize them into single-/multi-output diags)
void set_diagnostic(
    const std::vector<std::string> &all_diag_list,
    std::vector<std::string> &singleout_diag_list,
    std::vector<std::string> &multiout_diag_list,
    std::map<std::string, std::vector<std::string>> &diag_map);

// Create a diagnostic from a string representation of it.
// E.g., create the diag to compute fieldX_at_500hPa.
std::shared_ptr<AtmosphereDiagnostic>
create_diagnostic (const std::string& diag_name,
                   const std::shared_ptr<const AbstractGrid>& grid);
std::shared_ptr<AtmosphereDiagnostic>
create_diagnostic (const std::string& diag_name,
                   const std::shared_ptr<const AbstractGrid>& grid,
                   const std::vector<std::string>& multi_out_fields);

class AtmosDiagUtils {

public:
  std::map<std::string, std::regex> diag_regex;
  std::vector<std::string> multiout_diags;

  AtmosDiagUtils () {
    // Note: use grouping (the (..) syntax), so you can later query the content
    //       of each group in the matches output var!
    // Note: use raw string syntax R"(<string>)" to avoid having to escape the \ character
    // Note: the number for field_at_p/h can match positive integer/floating-point numbers
    diag_regex["FieldAtLevel"] = std::regex(R"(([A-Za-z0-9_]+)_at_(lev_(\d+)|model_(top|bot))$)");
    diag_regex["FieldAtPressureLevel"] = std::regex(R"(([A-Za-z0-9_]+)_at_(\d+(\.\d+)?)(hPa|mb|Pa)$)");
    diag_regex["FieldAtHeight"] = std::regex(R"(([A-Za-z0-9_]+)_at_(\d+(\.\d+)?)(m)_above_(sealevel|surface)$)");
    diag_regex["precip_surf_mass_flux"] = std::regex("precip_(liq|ice|total)_surf_mass_flux$");
    diag_regex["WaterPath"] = std::regex("(Ice|Liq|Rain|Rime|Vap)WaterPath$");
    diag_regex["NumberPath"] = std::regex("(Ice|Liq|Rain)NumberPath$");
    diag_regex["AeroComCld"] = std::regex("AeroComCld(Top|Bot)$");
    diag_regex["VaporFlux"] = std::regex("(Meridional|Zonal)VapFlux$");
    diag_regex["AtmBackTendDiag"] = std::regex("([A-Za-z0-9_]+)_atm_backtend$");
    diag_regex["PotentialTemperature"] = std::regex("(Liq)?PotentialTemperature$");
    diag_regex["VerticalLayer"] = std::regex("(z|geopotential|height)_(mid|int)$");
    diag_regex["HorizAvgDiag"] = std::regex("([A-Za-z0-9_]+)_horiz_avg$");

    // For now, we hardcore the multi-out diag candidates for ease
    multiout_diags.push_back("NumberPath");
  }
};

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
