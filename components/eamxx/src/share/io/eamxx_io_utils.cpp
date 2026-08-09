#include "share/io/eamxx_io_utils.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_config.hpp"

#include "share/io/eamxx_diag_spec.hpp"
#include "share/io/eamxx_legacy_diag_names.hpp"

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

std::shared_ptr<AbstractDiagnostic>
create_diagnostic (const std::string& diag_field_name,
                   const std::shared_ptr<const AbstractGrid>& grid)
{
  // A field computed by a diag that produces several of them is just a bare
  // name, with no expression structure to lower. Resolve it first, so that the
  // legacy translation below never gets a chance to read structure into it.
  const auto provider = DiagOutputsRegistry::instance().provider_of(diag_field_name);
  if (not provider.empty()) {
    ekat::ParameterList params(diag_field_name);
    params.set("grid_name",grid->name());
    return DiagnosticFactory::instance().create(provider,grid->get_comm(),params,grid);
  }

  // Translate a legacy mangled name (if that is what this is), parse it, and
  // lower it into the arguments the factory needs. See eamxx_diag_spec.hpp.
  // NOTE: this creates the diagnostic ALONE: its inputs are not set, and any
  //       sub-expression it depends on is not built. Use DiagBank to build a
  //       whole request, which is what the output streams do.
  const auto spec = lower_to_diag_spec(legacy_to_expr(diag_field_name),diag_field_name);

  auto params = spec.params;
  params.set("grid_name",grid->name());

  // An empty factory key means the name was not recognized as an expression, so
  // it can only be a diagnostic registered under that very name.
  const auto& key = spec.factory_key.empty() ? diag_field_name : spec.factory_key;

  return DiagnosticFactory::instance().create(key,grid->get_comm(),params,grid);
}

} // namespace scream
