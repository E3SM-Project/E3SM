#include "share/io/eamxx_io_utils.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_config.hpp"

#include <ekat_string_utils.hpp>

#include <fstream>
#include <regex>
#include <sstream>

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

util::TimeStamp parse_cf_time_units (const std::string& time_units)
{
  // Parse CF convention time units string (e.g., "seconds since 1970-01-01 00:00:00")
  // Expected format: "<unit> since <YYYY-MM-DD> <HH:MM:SS>"
  // where unit can be: seconds, minutes, hours, days
  
  // Find "since" keyword
  auto since_pos = time_units.find("since");
  EKAT_REQUIRE_MSG (since_pos != std::string::npos,
      "Error! Time units string does not contain 'since' keyword.\n"
      "  - time_units: " + time_units + "\n"
      "  - expected format: '<unit> since <YYYY-MM-DD> <HH:MM:SS>'\n");
  
  // Extract the datetime part after "since"
  std::string datetime_str = time_units.substr(since_pos + 5); // 5 = length of "since"
  
  // Trim leading whitespace
  size_t start = datetime_str.find_first_not_of(" \t\n\r");
  if (start != std::string::npos) {
    datetime_str = datetime_str.substr(start);
  }
  
  // Parse the datetime string
  // Expected format: "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" or variations
  int year = 0, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  
  // Find date and time parts
  auto space_pos = datetime_str.find(' ');
  std::string date_str = (space_pos != std::string::npos) ? datetime_str.substr(0, space_pos) : datetime_str;
  std::string time_str = (space_pos != std::string::npos) ? datetime_str.substr(space_pos + 1) : "";
  
  // Parse date part (YYYY-MM-DD)
  std::istringstream date_stream(date_str);
  char sep;
  date_stream >> year >> sep >> month >> sep >> day;
  
  EKAT_REQUIRE_MSG (date_stream.good() or date_stream.eof(),
      "Error! Failed to parse date from time units string.\n"
      "  - time_units: " + time_units + "\n"
      "  - date_str: " + date_str + "\n"
      "  - expected format: 'YYYY-MM-DD'\n");
  
  // Parse time part if present (HH:MM:SS)
  if (!time_str.empty()) {
    // Remove any trailing timezone information or extra characters
    auto tz_pos = time_str.find_first_not_of("0123456789:. ");
    if (tz_pos != std::string::npos) {
      time_str = time_str.substr(0, tz_pos);
    }
    
    // Trim trailing whitespace
    auto end = time_str.find_last_not_of(" \t\n\r");
    if (end != std::string::npos) {
      time_str = time_str.substr(0, end + 1);
    }
    
    if (!time_str.empty()) {
      std::istringstream time_stream(time_str);
      char sep1, sep2;
      time_stream >> hour;
      if (time_stream.peek() == ':') {
        time_stream >> sep1 >> minute;
        if (time_stream.peek() == ':') {
          time_stream >> sep2 >> second;
        }
      }
    }
  }
  
  // Create and return the timestamp
  return util::TimeStamp(year, month, day, hour, minute, second);
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
  std::regex field_at_l (R"()" + generic_field + R"(_at_(lev_(\d+)|model_(top|bot))$)");
  std::regex field_at_p (R"()" + generic_field + R"(_at_(\d+(\.\d+)?)(hPa|mb|Pa)$)");
  std::regex field_at_h (R"()" + generic_field + R"(_at_(\d+(\.\d+)?)(m)_above_(sealevel|surface)$)");
  std::regex surf_mass_flux ("precip_(liq|ice|total)_surf_mass_flux$");
  std::regex water_path ("(Ice|Liq|Rain|Rime|Vap)WaterPath$");
  std::regex number_path ("(Ice|Liq|Rain)NumberPath$");
  std::regex aerocom_cld ("AeroComCld(Top|Bot)$");
  std::regex vap_flux ("(Meridional|Zonal)VapFlux$");
  std::regex backtend (generic_field + "_atm_backtend$");
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
  } else if (std::regex_search(diag_field_name,matches,backtend)) {
    diag_name = "AtmBackTendDiag";
    params.set("grid_name",grid->name());
    params.set<std::string>("tendency_name",matches[1].str());
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
