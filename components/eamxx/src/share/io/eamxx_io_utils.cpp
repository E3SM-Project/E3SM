#include "share/io/eamxx_io_utils.hpp"

#include "share/io/eamxx_dexpr_diags.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_config.hpp"

#include <ekat_string_utils.hpp>

#include <fstream>
#include <map>
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

std::shared_ptr<AbstractDiagnostic>
create_diagnostic (const std::string& diag_field_name,
                   const std::shared_ptr<const AbstractGrid>& grid)
{
  // An expression, e.g. 'T_mid.interp(p_mid=500,units='hPa')' or '(qc+qv)*p_mid'.
  // A plain identifier comes back nullptr, and is resolved below.
  if (auto diag = dexpr_create_diagnostic(diag_field_name,grid)) {
    return diag;
  }

  // A plain identifier names a diag. A few diag classes serve a family of names,
  // where the name picks a variant (e.g., LiqWaterPath and IceWaterPath are both
  // WaterPath). Anything else is taken as the diag class name itself.
  struct Variant {
    std::string diag_name;
    std::string param_name;
    std::string param_value;
  };
  std::map<std::string,Variant> named_variants;
  for (std::string k : {"liq","ice","total"}) {
    named_variants["precip_" + k + "_surf_mass_flux"] = {"precip_surf_mass_flux","precip_type",k};
  }
  for (std::string k : {"Ice","Liq","Rain","Rime","Vap"}) {
    named_variants[k + "WaterPath"] = {"WaterPath","water_kind",k};
  }
  for (std::string k : {"Ice","Liq","Rain"}) {
    named_variants[k + "NumberPath"] = {"NumberPath","number_kind",k};
  }
  for (std::string k : {"Top","Bot"}) {
    named_variants["AeroComCld" + k] = {"AeroComCld","aero_com_cld_kind",k};
  }
  for (std::string k : {"Meridional","Zonal"}) {
    named_variants[k + "VapFlux"] = {"VaporFlux","wind_component",k};
  }
  named_variants["PotentialTemperature"]    = {"PotentialTemperature","temperature_kind","Tot"};
  named_variants["LiqPotentialTemperature"] = {"PotentialTemperature","temperature_kind","Liq"};

  std::string diag_name = diag_field_name;
  ekat::ParameterList params(diag_field_name);

  if (auto it=named_variants.find(diag_field_name); it!=named_variants.end()) {
    const auto& v = it->second;
    EKAT_REQUIRE_MSG (v.diag_name!="AeroComCld",
        "Error! AeroComCld diags are disabled for now. Contact developers.\n"
        "      Some recent development made the code produce bad values,\n"
        "      even runtime aborts due to NaNs.\n"
        "      An alternative is to request variables like cdnc_at_cldtop,\n"
        "      which remain unaffected and scientifically valid.\n");
    diag_name = v.diag_name;
    params.set<std::string>(v.param_name,v.param_value);
  } else if (diag_field_name=="dz") {
    diag_name = "VerticalLayer";
    params.set<std::string>("diag_name","dz");
    params.set<std::string>("vert_location","mid");
  } else {
    for (std::string n : {"z","geopotential","height"}) {
      for (std::string loc : {"mid","int"}) {
        if (diag_field_name==n + "_" + loc) {
          diag_name = "VerticalLayer";
          params.set<std::string>("diag_name",n);
          params.set<std::string>("vert_location",loc);
        }
      }
    }
  }

  auto comm = grid->get_comm();
  auto diag = DiagnosticFactory::instance().create(diag_name,comm,params,grid);

  return diag;
}

} // namespace scream
