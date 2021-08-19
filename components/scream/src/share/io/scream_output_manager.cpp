#include "scream_output_manager.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<const fm_type>& field_mgr,
       const bool runtype_restart)
{
  m_io_comm   = io_comm;
  m_params    = params;
  m_field_mgr = field_mgr;
  m_runtype_restart = runtype_restart;

  // Check for restart output
  // Restarts are a special kind of output that have unique characteristics.
  // The following constructs the parameter list to control the restart output stream.
  // The user can trigger restart output at runtime by adding the sublist
  // "Restart Control" to the SCREAM control YAML file.
  // A typical restart control entry may look like:
  // ----
  // Restart Contol
  //   FREQUENCY
  //     OUT_N: INT
  //     OUT_OPTION: STRING
  // ----
  // where OUT_OPTION is the units of the output frequency (ex. Steps, Months, etc.),
  // and OUT_N>0 is the actual frequency. E.g., OUT_N=2 meanse every other ${OUT_OPTION}

  // The output history should have the same restart frequency as the model restart output.
  int restart_freq = 0;
  std::string restart_freq_units = "None";
  if (m_params.isSublist("Restart Control")) {
    // Get restart parameters, and create a param list for the model restart output
    auto& out_params = m_params.sublist("Restart Control");
    make_restart_param_list(out_params);
    new_output(out_params,true);

    // Gather restart frequency info for restart history files
    auto& freq_params = out_params.sublist("FREQUENCY");
    restart_freq = freq_params.get<Int>("OUT_N");
    restart_freq_units = freq_params.get<std::string>("OUT_OPTION");
  }

  // Construct and store an output stream instance for each output request.
  // Typical output is controlled by parameter list which is stored in individual YAML files.
  // A list of those files is passed to the output manager via the parameter list using the
  // key "Output YAML Files".  See srd/share/io/tests for examples.
  std::vector<std::string> empty;
  auto& list_of_files = m_params.get("Output YAML Files",empty);
  for (auto& it : list_of_files) {
    ekat::ParameterList out_params(it);
    parse_yaml_file(it,out_params);
    out_params.set("CHECKPOINT FREQUENCY",restart_freq);
    out_params.set("OUT_OPTION",restart_freq_units);
    new_output(out_params,false);
  }
}
/*===============================================================================================*/
void OutputManager::run(util::TimeStamp& current_ts)
{
  for (auto& it : m_output_streams) {
    it->run(current_ts);
  }
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  for (auto& it : m_output_streams) {
    it->finalize();
  }
  m_output_streams.clear();
  m_io_comm = ekat::Comm();
  m_params = ekat::ParameterList("");
  m_field_mgr = nullptr;
  m_runtype_restart = false;
}
/*===============================================================================================*/
void OutputManager::
new_output(const ekat::ParameterList& params, const bool model_restart_output)
{
  auto output = std::make_shared<output_type>(m_io_comm,params,m_field_mgr,
                                              m_runtype_restart, model_restart_output);
  m_output_streams.push_back(output);
}
/*===============================================================================================*/
void OutputManager::make_restart_param_list(ekat::ParameterList& params)
{
  // Given the unique nature of restart files, this routine sets up the specific requirements.

  //TODO change this to some sort of case_name, TODO: control so suffix is .r.nc
  params.set<std::string>("FILENAME", "scorpio_restart_test");

  // Parse the parameters that controls this output instance.
  // Model restart themselves are instant snapshots (no averaging of any kind).
  params.set<std::string>("AVERAGING TYPE","Instant");

  if (not params.isParameter("GRID")) {
    params.set<std::string>("GRID",m_field_mgr->get_grid()->name());
  }
  auto& freq_params = params.sublist("FREQUENCY");
  freq_params.get<Int>("OUT_MAX_STEPS",1);

  // Grab the restart fields from the field manager.
  // If a developer wants a field to be stored in the restart file than they must add the "restart" group tag
  // to that field when registering the field with the field manager.
  auto restart_fields = m_field_mgr->get_field_group("restart");

  // WARNING: the vector restart_fields.m_info->m_fields_names contains CaseInsensitiveString's.
  //          This can be a problem, since virtually all places will try to extract a
  //          std::vector<std::string> from the parameter list, causing a any_cast error.
  //          To fix this, create a std::vector<std::string> on the fly
  std::vector<std::string> fields_names;
  for (const auto& fn : restart_fields.m_info->m_fields_names) {
    fields_names.push_back(fn);
  }
  params.set("FIELDS",fields_names);
}
/*===============================================================================================*/
} // namespace scream
