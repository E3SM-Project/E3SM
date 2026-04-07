#include "share/io/eamxx_io_utils.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_config.hpp"

#include <ekat_string_utils.hpp>

#include <fstream>
#include <regex>

namespace scream {

std::string find_filename_in_rpointer (
    const std::string& filename_prefix,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0,
    const bool allow_not_found,
    const OutputAvgType avg_type,
    const IOControl& control)
{
  std::string filename;
  bool found = false;
  std::string content;
  std::string suffix = model_restart ? ".r." : ".rhist.";
  std::string pattern_str = filename_prefix + suffix;

  // The AD will pass a default constructed control, since it doesn't know the values
  // of REST_N/REST_OPTION used in the previous run. Also, model restart is *always* INSTANT.
  if (model_restart) {
    EKAT_REQUIRE_MSG (avg_type==OutputAvgType::Instant,
        "Error! Model restart output should have INSTANT avg type.\n"
        " - input avg_type: " + e2str(avg_type) + "\n");
    pattern_str += e2str(OutputAvgType::Instant) + R"(.n(step|sec|min|hour|day|month|year)s_x\d+)";
  } else {
    EKAT_REQUIRE_MSG (control.output_enabled(),
        "Error! When restarting an output stream, we need a valid IOControl structure.\n"
        " - filename prefix: " + filename_prefix + "\n");
    pattern_str += e2str(avg_type) + "." + control.frequency_units + "_x" + std::to_string(control.frequency);
  }
  if (is_scream_standalone()) {
    pattern_str += ".np" + std::to_string(comm.size());
  }
  pattern_str += "." + run_t0.to_string() + ".nc";
  std::regex pattern (pattern_str);

  if (comm.am_i_root()) {
    std::ifstream rpointer_file;

    std::string line;
    rpointer_file.open("rpointer.atm");

    while (std::getline(rpointer_file,line)) {
      content += line + "\n";

      if (std::regex_search(line,pattern)) {
        filename = line;
        found = true;
        break;
      }
    }
  }

  int ifound = int(found);
  comm.broadcast(&ifound,1,0);
  found = bool(ifound);

  if (found) {
    // Have the root rank communicate the nc filename
    broadcast_string(filename,comm,comm.root_rank());
  } else if (not allow_not_found) {
    broadcast_string(content,comm,comm.root_rank());

    if (model_restart) {
      EKAT_ERROR_MSG (
          "Error! Restart requested, but no model restart file found in 'rpointer.atm'.\n"
          "   model restart filename prefix: " + filename_prefix + "\n"
          "   model restart filename pattern: " + pattern_str + "\n"
          "   run t0           : " + run_t0.to_string() + "\n"
          "   rpointer content:\n" + content + "\n\n");
    } else {
      EKAT_ERROR_MSG (
          "Error! Restart requested, but no history restart file found in 'rpointer.atm'.\n"
          "   hist restart filename prefix: " + filename_prefix + "\n"
          "   hist restart filename pattern: " + pattern_str + "\n"
          "   run t0           : " + run_t0.to_string() + "\n"
          "   avg_type         : " + e2str(avg_type) + "\n"
          "   output freq      : " + std::to_string(control.frequency) + "\n"
          "   output freq units: " + control.frequency_units + "\n"
          "   rpointer content:\n" + content + "\n\n"
          " Did you change output specs (avg type, freq, or freq units) across restart? If so, please, remember that it is not allowed.\n"
          " It is also possible you are using a rhist file create before commit 6b7d441330d. That commit changed how rhist file names\n"
          " are formed. In particular, we no longer use INSTANT.${REST_OPTION}_x${REST_N}, but we use the avg type, and freq/freq_option\n"
          " of the output stream (to avoid name clashes if 2 streams only differ for one of those). If you want to use your rhist file,\n"
          " please rename it, so that the avg-type, freq, and freq_option reflect those of the output stream.\n");
    }
  }

  return filename;
}

void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps)
{
  scorpio::set_attribute(filename,"GLOBAL",ts_name,ts.to_string());
  if (write_nsteps) {
    scorpio::set_attribute(filename,"GLOBAL",ts_name+"_nsteps",ts.get_num_steps());
  }
}

util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps)
{
  auto ts = util::str_to_time_stamp(scorpio::get_attribute<std::string>(filename,"GLOBAL",ts_name));
  if (read_nsteps and scorpio::has_attribute(filename,"GLOBAL",ts_name+"_nsteps")) {
    ts.set_num_steps(scorpio::get_attribute<int>(filename,"GLOBAL",ts_name+"_nsteps"));
  }
  return ts;
}

std::pair<util::TimeStamp,int>
parse_cf_time_units (const std::string& units_str)
{
  // Find the " since " separator
  const std::string sep = " since ";
  auto pos = units_str.find(sep);
  EKAT_REQUIRE_MSG (pos!=std::string::npos,
      "Error! Could not parse CF-compliant time units string.\n"
      " - units string: '" + units_str + "'\n"
      " - Expected format: '<unit> since <YYYY-MM-DD> [HH[:MM[:SS]]]'\n");

  // Map unit name to seconds using exact matching
  auto unit_str = units_str.substr(0,pos);
  int time_mult;
  if (unit_str == "second" || unit_str == "seconds") {
    time_mult = 1;
  } else if (unit_str == "minute" || unit_str == "minutes") {
    time_mult = 60;
  } else if (unit_str == "hour" || unit_str == "hours") {
    time_mult = 3600;
  } else if (unit_str == "day" || unit_str == "days") {
    time_mult = 86400;
  } else {
    EKAT_ERROR_MSG (
        "Error! Unsupported time unit in CF-compliant units string.\n"
        " - units string: '" + units_str + "'\n"
        " - unit: '" + unit_str + "'\n"
        " - Supported units: second(s), minute(s), hour(s), day(s)\n");
  }

  // Extract and parse the reference date/time string
  auto date_str = units_str.substr(pos + sep.size());

  // First try EAMxx internal format: YYYY-MM-DD-SSSSS
  auto ref_ts = util::str_to_time_stamp(date_str);
  if (ref_ts.is_valid()) {
    return {ref_ts, time_mult};
  }

  // Try CF format: YYYY-MM-DD [HH[:MM[:SS]]]
  EKAT_REQUIRE_MSG (date_str.size()>=10 && date_str[4]=='-' && date_str[7]=='-',
      "Error! Could not parse reference date in CF-compliant time units string.\n"
      " - units string: '" + units_str + "'\n"
      " - date/time part: '" + date_str + "'\n");

  int yy = 0, mon = 0, dd = 0;
  try {
    yy  = std::stoi(date_str.substr(0,4));
    mon = std::stoi(date_str.substr(5,2));
    dd  = std::stoi(date_str.substr(8,2));
  } catch (...) {
    EKAT_ERROR_MSG (
        "Error! Could not parse date components in CF-compliant time units string.\n"
        " - units string: '" + units_str + "'\n"
        " - date/time part: '" + date_str + "'\n");
  }
  int hh = 0, min = 0, sec = 0;

  if (date_str.size() > 10) {
    // Expect a space or 'T' separator, then HH[:MM[:SS]]
    char sep_char = date_str[10];
    EKAT_REQUIRE_MSG (sep_char==' ' || sep_char=='T',
        "Error! Invalid separator between date and time in CF-compliant units string.\n"
        " - units string: '" + units_str + "'\n"
        " - date/time part: '" + date_str + "'\n");
    auto time_part = date_str.substr(11);
    EKAT_REQUIRE_MSG (time_part.size() >= 2,
        "Error! Could not parse time in CF-compliant units string.\n"
        " - units string: '" + units_str + "'\n");
    try {
      hh = std::stoi(time_part.substr(0,2));
      if (time_part.size() >= 5 && time_part[2] == ':') {
        min = std::stoi(time_part.substr(3,2));
        if (time_part.size() >= 8 && time_part[5] == ':') {
          sec = std::stoi(time_part.substr(6,2));
        }
      }
    } catch (...) {
      EKAT_ERROR_MSG (
          "Error! Could not parse time components in CF-compliant time units string.\n"
          " - units string: '" + units_str + "'\n"
          " - time part: '" + time_part + "'\n");
    }
  }

  ref_ts = util::TimeStamp({yy,mon,dd},{hh,min,sec});
  EKAT_REQUIRE_MSG (ref_ts.is_valid(),
      "Error! Could not create a valid TimeStamp from CF-compliant units string.\n"
      " - units string: '" + units_str + "'\n");

  return {ref_ts, time_mult};
}

std::shared_ptr<AtmosphereDiagnostic>
create_diagnostic (const std::string& diag_field_name,
                   const std::shared_ptr<const AbstractGrid>& grid)
{
  // Note: use grouping (the (..) syntax), so you can later query the content
  //       of each group in the matches output var!
  // Note: use raw string syntax R"(<string>)" to avoid having to escape the \ character
  // Note: the number for field_at_p/h can match positive integer/floating-point numbers
  // Start with a generic for a field name allowing for all letters, all numbers, dash, dot, plus, minus, product, and division
  // Escaping all the special ones just in case
  std::string generic_field = "([A-Za-z0-9_.+\\-\\*\\÷]+)";

  // ── Built-in aliases ──────────────────────────────────────────────────────
  // Recognized shorthand patterns expand to canonical composable expressions.
  // Checked before all other regexes; expansion recurses into create_diagnostic.
  {
    std::smatch alias_matches;
    std::regex bt (generic_field + "_atm_backtend$");
    if (std::regex_search(diag_field_name, alias_matches, bt)) {
      const auto f = alias_matches[1].str();
      return create_diagnostic(f + "_minus_" + f + "_prev_over_dt", grid);
    }
    // More built-in aliases may be added here.
  }
  std::regex field_at_l (R"()" + generic_field + R"(_at_(lev_(\d+)|model_(top|bot))$)");
  std::regex field_at_p (R"()" + generic_field + R"(_at_(\d+(\.\d+)?)(hPa|mb|Pa)$)");
  std::regex field_at_h (R"()" + generic_field + R"(_at_(\d+(\.\d+)?)(m)_above_(sealevel|surface)$)");
  std::regex surf_mass_flux ("precip_(liq|ice|total)_surf_mass_flux$");
  std::regex water_path ("(Ice|Liq|Rain|Rime|Vap)WaterPath$");
  std::regex number_path ("(Ice|Liq|Rain)NumberPath$");
  std::regex aerocom_cld ("AeroComCld(Top|Bot)$");
  std::regex vap_flux ("(Meridional|Zonal)VapFlux$");
  std::regex field_prev (generic_field + "_prev$");
  std::regex field_over_dt (generic_field + "_over_dt$");
  std::regex pot_temp ("(Liq)?PotentialTemperature$");
  std::regex vert_layer ("(z|geopotential|height)_(mid|int)$");
  std::regex horiz_avg (generic_field + "_horiz_avg$");
  std::regex vert_contract (generic_field + "_vert_(avg|sum)(_((dp|dz)_weighted))?$");
  std::regex zonal_avg (R"()" + generic_field + R"(_zonal_avg_(\d+)_bins$)");
  std::regex conditional_sampling (R"()" + generic_field + R"(_where_)" + generic_field + R"(_(gt|ge|eq|ne|le|lt)_([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)$)");
  std::regex binary_ops (generic_field + "_" "(plus|minus|times|over)" + "_" + generic_field + "$");
  std::regex histogram (R"()" + generic_field + R"(_histogram_(\d+(\.\d+)?(_\d+(\.\d+)?)+)$)");
  std::regex vert_derivative (generic_field + "_(p|z)vert_derivative$");

  std::string diag_name;
  std::smatch matches;
  ekat::ParameterList params(diag_field_name);

  if (std::regex_search(diag_field_name,matches,field_at_l)) {
    params.set("field_name",matches[1].str());
    params.set("grid_name",grid->name());
    params.set("vertical_location", matches[2].str());
    diag_name = "FieldAtLevel";
  } else if (std::regex_search(diag_field_name,matches,field_at_p)) {
    params.set("field_name",matches[1].str());
    params.set("grid_name",grid->name());
    params.set("pressure_value",matches[2].str());
    params.set("pressure_units", matches[4].str());
    diag_name = "FieldAtPressureLevel";
  } else if (std::regex_search(diag_field_name,matches,field_at_h)) {
    params.set("field_name",matches[1].str());
    params.set("grid_name",grid->name());
    params.set("height_value",matches[2].str());
    params.set("height_units",matches[4].str());
    params.set("surface_reference", matches[5].str());
    diag_name = "FieldAtHeight";
  } else if (std::regex_search(diag_field_name,matches,surf_mass_flux)) {
    diag_name = "precip_surf_mass_flux";
    params.set<std::string>("precip_type",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,water_path)) {
    diag_name = "WaterPath";
    params.set<std::string>("water_kind",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,number_path)) {
    diag_name = "NumberPath";
    params.set<std::string>("number_kind",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,aerocom_cld)) {
    EKAT_ERROR_MSG("Error! AeroComCld diags are disabled for now. Contact developers.\n"
                    "      Some recent development made the code produce bad values,\n"
                    "      even runtime aborts due to NaNs.\n"
                    "      An alternative is to request variables like cdnc_at_cldtop,\n"
                    "      which remain unaffected and scientifically valid.\n");
    diag_name = "AeroComCld";
    params.set<std::string>("aero_com_cld_kind",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,vap_flux)) {
    diag_name = "VaporFlux";
    params.set<std::string>("wind_component",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,field_over_dt)) {
    // NOTE: _over_dt must be checked before binary_ops to prevent "X_over_dt" from being
    // misinterpreted as BinaryOpsDiag(X, over, dt).
    diag_name = "FieldOverDtDiag";
    params.set("grid_name",grid->name());
    params.set<std::string>("field_name",matches[1].str());
  } else if (std::regex_search(diag_field_name,matches,pot_temp)) {
    diag_name = "PotentialTemperature";
    params.set<std::string>("temperature_kind", matches[1].str()!="" ? matches[1].str() : std::string("Tot"));
  } else if (std::regex_search(diag_field_name,matches,vert_layer)) {
    diag_name = "VerticalLayer";
    params.set<std::string>("diag_name",matches[1].str());
    params.set<std::string>("vert_location",matches[2].str());
  } else if (diag_field_name=="dz") {
    diag_name = "VerticalLayer";
    params.set<std::string>("diag_name","dz");
    params.set<std::string>("vert_location","mid");
  }
  else if (std::regex_search(diag_field_name,matches,horiz_avg)) {
    diag_name = "HorizAvgDiag";
    params.set("grid_name",grid->name());
    params.set<std::string>("field_name",matches[1].str());
  }
  else if (std::regex_search(diag_field_name,matches,vert_contract)) {
    diag_name = "VertContractDiag";
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", matches[1].str());
    params.set<std::string>("contract_method", matches[2].str());
    // The 3rd match an optional _(dp|dz)_weighted, so check if it was matched
    if (matches[3].matched) {
      // note that the 4th match is (dp|dz)_weighted, while the 5th is (dp|dz)
      params.set<std::string>("weighting_method", matches[5].str());
    }
  }
  else if (std::regex_search(diag_field_name,matches,vert_derivative)) {
    diag_name = "VertDerivativeDiag";
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", matches[1].str());
    params.set<std::string>("derivative_method", matches[2].str());
  }
  else if (std::regex_search(diag_field_name,matches,zonal_avg)) {
    diag_name = "ZonalAvgDiag";
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", matches[1].str());
    params.set<std::string>("number_of_zonal_bins", matches[2].str());
  }
  else if (std::regex_search(diag_field_name,matches,conditional_sampling)) {
    diag_name = "ConditionalSampling";
    params.set("grid_name", grid->name());
    params.set<std::string>("input_field", matches[1].str());
    params.set<std::string>("condition_field", matches[2].str());
    params.set<std::string>("condition_operator", matches[3].str());
    params.set<std::string>("condition_value", matches[4].str());
  }
  else if (std::regex_search(diag_field_name,matches,binary_ops)) {
    diag_name = "BinaryOpsDiag";
    params.set("grid_name", grid->name());
    params.set<std::string>("arg1", matches[1].str());
    params.set<std::string>("arg2", matches[3].str());
    params.set<std::string>("binary_op", matches[2].str());
  }
  // NOTE: field_prev must be checked AFTER binary_ops so that "X_minus_X_prev"
  // is parsed as BinaryOpsDiag(X, minus, X_prev) rather than FieldPrevDiag(X_minus_X).
  else if (std::regex_search(diag_field_name,matches,field_prev)) {
    diag_name = "FieldPrevDiag";
    params.set("grid_name",grid->name());
    params.set<std::string>("field_name",matches[1].str());
  }
  else if (std::regex_search(diag_field_name,matches,histogram)) {
    diag_name = "HistogramDiag";
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", matches[1].str());
    params.set<std::string>("bin_configuration", matches[2].str());
  }
  else
  {
    // No existing special regex matches, so we assume that the diag field name IS the diag name.
    diag_name = diag_field_name;
  }

  auto comm = grid->get_comm();
  auto diag = AtmosphereDiagnosticFactory::instance().create(diag_name,comm,params);
  auto gm = std::make_shared<LibraryGridsManager>(grid);
  diag->set_grids(gm);

  return diag;
}

} // namespace scream
