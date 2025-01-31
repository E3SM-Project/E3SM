#include "share/io/scream_io_utils.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/library_grids_manager.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_config.hpp"

#include <fstream>

namespace scream {

std::string find_filename_in_rpointer (
    const std::string& filename_prefix,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0,
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

  if (not found) {
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

  // Have the root rank communicate the nc filename
  broadcast_string(filename,comm,comm.root_rank());

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

void set_diagnostic(const std::vector<std::string> &all_diag_list,
                    std::vector<std::string> &singleout_diag_list,
                    std::vector<std::string> &multiout_diag_list,
                    std::map<std::string, std::vector<std::string>> &diag_map) {
  AtmosDiagUtils diag_utils;
  auto diag_regex     = diag_utils.diag_regex;
  auto multiout_diags = diag_utils.multiout_diags;

  for(const auto &name : all_diag_list) {
    bool is_multi = false;
    for(const auto &diag : multiout_diags) {
      if(std::regex_match(name, diag_regex[diag])) {
        if(std::find(multiout_diag_list.begin(), multiout_diag_list.end(),
                     diag) == multiout_diag_list.end()) {
          multiout_diag_list.push_back(diag);
        }
        diag_map[diag].push_back(name);
        is_multi = true;
        break;
      }
    }
    if(!is_multi) {
      singleout_diag_list.push_back(name);
      // empty for single-output diagnostic calls
      diag_map[name] = std::vector<std::string>();
    }
  }
}

std::shared_ptr<AtmosphereDiagnostic> create_diagnostic(
    const std::string &diag_field_name,
    const std::shared_ptr<const AbstractGrid> &grid) {
  return create_diagnostic(diag_field_name, grid, std::vector<std::string>());
}

std::shared_ptr<AtmosphereDiagnostic> create_diagnostic(
    const std::string &diag_field_name,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::vector<std::string> &multi_out_fields) {
  AtmosDiagUtils diag_utils;
  auto diag_regex     = diag_utils.diag_regex;
  auto multiout_diags = diag_utils.multiout_diags;

  std::string diag_name;
  std::smatch matches;
  ekat::ParameterList params(diag_field_name);

  // Take care of the multi-out diags first

  if(std::regex_search(diag_field_name, matches, diag_regex["FieldAtLevel"])) {
    params.set("field_name", matches[1].str());
    params.set("grid_name", grid->name());
    params.set("vertical_location", matches[2].str());
    diag_name = "FieldAtLevel";
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["FieldAtPressureLevel"])) {
    params.set("field_name", matches[1].str());
    params.set("grid_name", grid->name());
    params.set("pressure_value", matches[2].str());
    params.set("pressure_units", matches[4].str());
    diag_name = "FieldAtPressureLevel";
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["FieldAtHeight"])) {
    params.set("field_name", matches[1].str());
    params.set("grid_name", grid->name());
    params.set("height_value", matches[2].str());
    params.set("height_units", matches[4].str());
    params.set("surface_reference", matches[5].str());
    diag_name = "FieldAtHeight";
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["precip_surf_mass_flux"])) {
    diag_name = "precip_surf_mass_flux";
    params.set<std::string>("precip_type", matches[1].str());
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["WaterPath"])) {
    diag_name = "WaterPath";
    params.set<std::string>("Water Kind", matches[1].str());
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["NumberPath"])) {
    // TODO: implement smarter way to create diags with multiple outputs
    if(multi_out_fields.size() == 0) {
      // This is not really used, but leave it just in case
      diag_name = "NumberPath";
      params.set<std::vector<std::string>>(
          "Number Kinds", std::vector<std::string>{matches[1].str()});
    } else {
      for(const auto &multi_out_field : multi_out_fields) {
        std::vector<std::string> multi_out_inputs;
        if(std::regex_search(multi_out_field, matches,
                             diag_regex["NumberPath"])) {
          diag_name = "NumberPath";
          multi_out_inputs.push_back(matches[1].str());
        }
      }
      params.set<std::vector<std::string>>("Number Kinds", multi_out_fields);
    }
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["AeroComCld"])) {
    diag_name = "AeroComCld";
    params.set<std::string>("AeroComCld Kind", matches[1].str());
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["VaporFlux"])) {
    diag_name = "VaporFlux";
    params.set<std::string>("Wind Component", matches[1].str());
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["AtmBackTendDiag"])) {
    diag_name = "AtmBackTendDiag";
    // Set the grid_name
    params.set("grid_name", grid->name());
    params.set<std::string>("Tendency Name", matches[1].str());
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["PotentialTemperature"])) {
    diag_name = "PotentialTemperature";
    params.set<std::string>("Temperature Kind", matches[1].str() != ""
                                                    ? matches[1].str()
                                                    : std::string("Tot"));
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["VerticalLayer"])) {
    diag_name = "VerticalLayer";
    params.set<std::string>("diag_name", matches[1].str());
    params.set<std::string>("vert_location", matches[2].str());
  } else if(diag_field_name == "dz") {
    diag_name = "VerticalLayer";
    params.set<std::string>("diag_name", "dz");
    params.set<std::string>("vert_location", "mid");
  } else if(std::regex_search(diag_field_name, matches,
                              diag_regex["HorizAvgDiag"])) {
    diag_name = "HorizAvgDiag";
    // Set the grid_name
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", matches[1].str());
  } else {
    // No existing special regex matches, so we assume that the diag field name
    // IS the diag name.
    diag_name = diag_field_name;
  }

  auto comm = grid->get_comm();
  auto diag =
      AtmosphereDiagnosticFactory::instance().create(diag_name, comm, params);
  auto gm = std::make_shared<LibraryGridsManager>(grid);
  diag->set_grids(gm);

  return diag;
}

}  // namespace scream
