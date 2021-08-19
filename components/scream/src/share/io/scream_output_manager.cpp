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

  // Check for model restart output
  // Restarts are a special kind of output that have unique characteristics.
  // The following constructs the parameter list to control the restart output stream.
  // The user can trigger restart output at runtime by adding the sublist
  // "Model Restart" to the SCREAM control YAML file.
  // A typical restart control entry may look like:
  // ----
  // Model Restart:
  //   Output:
  //     Frequency:       INT
  //     Frequency Units: STRING
  // ----
  // where Frequency Units is the units of the output frequency (ex. Steps, Months, etc.),
  // and Frequency>0 is the actual frequency. E.g., Frequency=2 meanse every other ${Frequency Units}

  // The output history should have the same restart frequency as the model restart output.
  int checkpoint_freq = 0;
  std::string checkpoint_freq_units = "None";
  if (m_params.isSublist("Model Restart")) {
    // Get restart parameters, and create a param list for the model restart output
    auto& out_params = m_params.sublist("Model Restart");
    make_restart_param_list(out_params);
    add_output_stream(out_params,true);

    // Gather restart frequency info, since the checkpointing info must match them
    checkpoint_freq       = out_params.sublist("Output").get<Int>("Frequency");
    checkpoint_freq_units = out_params.sublist("Output").get<std::string>("Frequency Units");
  }

  // Construct and store an output stream instance for each output request.
  // Typical output is controlled by parameter list which is stored in individual YAML files.
  // A list of those files is passed to the output manager via the parameter list.
  // In particular, check the 'Output YAML Files' sublist, and look for the entry
  // with the same name as the grid.
  // See srd/share/io/tests for examples.
  std::vector<std::string> empty;
  const auto& grid_name = m_field_mgr->get_grid()->name();
  auto& list_of_files = m_params.sublist("Output YAML Files").get(grid_name,empty);
  for (auto& it : list_of_files) {
    ekat::ParameterList out_params(it);
    parse_yaml_file(it,out_params);
    auto& checkpointing_params = out_params.sublist("Checkpointing");
    EKAT_REQUIRE_MSG (checkpointing_params.get("Frequency",checkpoint_freq)==checkpoint_freq,
        "Error! Output checkpointing frequency is different from the model restart output one.\n"
        "       In case model restart output is generated, checkpointing must have the same frequency.\n");
    EKAT_REQUIRE_MSG (checkpointing_params.get("Frequency Units",checkpoint_freq_units)==checkpoint_freq_units,
        "Error! Output checkpointing frequency units are different from the model restart output one.\n"
        "       In case model restart output is generated, checkpointing must have the same units.\n");
    add_output_stream(out_params,false);
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
add_output_stream(const ekat::ParameterList& params, const bool model_restart_output)
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
  params.set<std::string>("Casename", "scorpio_restart_test");

  // Parse the parameters that controls this output instance.
  // Model restart themselves are instant snapshots (no averaging of any kind).
  params.set<std::string>("Averaging Type","Instant");

  params.set<Int>("Max Snapshots Per File",1);

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
  params.set("Fields",fields_names);
}
/*===============================================================================================*/
} // namespace scream
