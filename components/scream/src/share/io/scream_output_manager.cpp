#include "scream_output_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <fstream>
#include <memory>

namespace scream
{

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& t0,
       const bool is_restarted_run)
{
  using map_t = std::map<std::string,std::shared_ptr<fm_type>>;
  map_t fms;
  fms[field_mgr->get_grid()->name()] = field_mgr;
  setup(io_comm,params,fms,grids_mgr,t0,is_restarted_run);
}

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& t0,
       const bool is_restarted_run)
{
  m_io_comm = io_comm;
  m_t0      = t0;
  m_is_restarted_run = is_restarted_run;

  // Check for model restart output
  set_params(params,field_mgrs);

  // Output control
  auto& out_control_pl = m_params.sublist("Output Control");
  m_output_control.frequency  = out_control_pl.get<int>("Frequency");
  m_output_control.frequency_units = out_control_pl.get<std::string>("Frequency Units");
  m_output_control.nsteps_since_last_write = 0;

  // File specs
  m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File");
  m_output_file_specs.num_snapshots_in_file = 0;
  m_output_file_specs.filename_with_time_string = out_control_pl.get("Timestamp in Filename",true);
  m_output_file_specs.filename_with_mpiranks    = out_control_pl.get("MPI Ranks in Filename",false);

  // For each grid, create a separate output stream.
  const auto& fields_pl = m_params.sublist("Fields");
  for (auto it=fields_pl.params_names_cbegin(); it!=fields_pl.params_names_cend(); ++it) {
    const auto& gname = *it;
    EKAT_REQUIRE_MSG (grids_mgr->has_grid(gname),
        "Error! Output requested on grid '" + gname + "', but the grids manager does not store such grid.\n");

    EKAT_REQUIRE_MSG (field_mgrs.find(gname)!=field_mgrs.end(),
        "Error! Output requested on grid '" + gname + "', but no field manager is available for such grid.\n");

    auto stream_params = params;
    stream_params.set("Grid",gname);
    auto output = std::make_shared<output_type>(m_io_comm,stream_params,field_mgrs.at(gname),grids_mgr);
    m_output_streams.push_back(output);
  }

  // Check if we need to restart the output history
  const auto has_restart_data = (m_avg_type!="INSTANT");
  if (has_restart_data) {
    // This avg_type needs to save some info in order to restart the output.
    // E.g., we might save 30-day avg value for field F, but due to job size
    // break the run into three 10-day runs. We then need to save the state of
    // our averaging in a "restart" file (e.g., the current avg and avg_count).
    // Note: the user might decide *not* to restart the output, so give the option
    //       of disabling the restart. Also, the user might want to change the
    //       casename, so allow to specify a different casename for the restart file.
    bool perform_restart = is_restarted_run;
    std::string hist_restart_casename;

    if (m_params.isSublist("Restart")) {
      const auto& pl = params.sublist("Restart");
      perform_restart &= pl.isParameter("Perform Restart") ?
                         pl.get<bool>("Perform Restart") : true;;
      hist_restart_casename = pl.isParameter("Casename") ? 
                              pl.get<std::string>("Casename") : m_casename;
    }

    if (m_params.isSublist("Checkpoint Control")) {
      // Output control
      auto& pl = m_params.sublist("Checkpoint Control"); 
      m_checkpoint_control.frequency  = pl.get<int>("Frequency");
      m_checkpoint_control.frequency_units = pl.get<std::string>("Frequency Units");
      m_checkpoint_control.nsteps_since_last_write = 0;

      // File specs
      m_checkpoint_file_specs.max_snapshots_in_file = 1;
      m_checkpoint_file_specs.num_snapshots_in_file = 0;
      m_checkpoint_file_specs.filename_with_time_string = pl.get("Timestamp in Filename",true);
      m_checkpoint_file_specs.filename_with_mpiranks    = pl.get("MPI Ranks in Filename",false);
    }

    if (perform_restart) {
      bool found = false;
      std::string filename;
      if (m_io_comm.am_i_root()) {
        std::ifstream rpointer_file;
        rpointer_file.open("rpointer.atm");
        while (rpointer_file >> filename) {
          if (filename.find(hist_restart_casename) != std::string::npos &&
              filename.find("rhist") != std::string::npos) {
            found = true;
            break;
          }
        }
      }
      m_io_comm.broadcast(&found,1,0);

      // If the history restart file is not found, it must be because the last
      // model restart step coincided with a model output step, in which case
      // a restart history file is not written.
      // If that's the case, *disable* output restart, by setting
      //   'Restart'->'Perform Restart' = false
      // in the input parameter list
      EKAT_REQUIRE_MSG (found,
          "Error! Output restart requested, but the no history restart file found in 'rpointer.atm'.\n"
          "   restart file name root: " + hist_restart_casename + "\n");

      // Have the root rank communicate the nc filename
      broadcast_string(filename,m_io_comm,m_io_comm.root_rank());

      // Restart each stream
      for (auto stream : m_output_streams) {
        stream->restart(filename);
      }

      // Restart the output control
      ekat::ParameterList res_params("Input Parameters");
      res_params.set<std::string>("Filename",filename);
      AtmosphereInput hist_restart (m_io_comm,res_params);
      m_output_control.nsteps_since_last_write = hist_restart.read_int_scalar("avg_count");
    }
  }
}
/*===============================================================================================*/
void OutputManager::run(util::TimeStamp& timestamp)
{
  using namespace scorpio;

  // Check if we need to open a new file
  ++m_output_control.nsteps_since_last_write;
  ++m_checkpoint_control.nsteps_since_last_write;

  // Check if this is a write step (and what kind)
  const bool is_output_step     = m_output_control.is_write_step();
  const bool is_checkpoint_step = m_checkpoint_control.is_write_step() && not is_output_step;
  const bool is_write_step      = is_output_step || is_checkpoint_step;

  // If neither output or checkpoint, these won't be used anyways,
  // so no need to check if is_write_step == true.
  auto& control   = is_checkpoint_step ? m_checkpoint_control : m_output_control;
  auto& filespecs = is_checkpoint_step ? m_checkpoint_file_specs : m_output_file_specs;
  auto& filename  = filespecs.filename;

  // Compute filename (if write step)
  if (is_write_step) {
    // Check if we need to open a new file
    if (not filespecs.is_open) {
      // Compute new file name
      filename = compute_filename_root(control);
      if (m_params.get("Timestamp in Filename", true)) {
        filename += "." + timestamp.to_string();
      }
      if (is_output_step) {
        filename += m_is_model_restart_output ? ".r.nc" : ".nc";
      } else if (is_checkpoint_step) {
        filename += ".rhist.nc";
      } else {
        filename += ".nc";
      }

      // Register new netCDF file for output. First, check no other output managers
      // are trying to write on the same file
      EKAT_REQUIRE_MSG (not is_file_open_c2f(filename.c_str(),Write),
          "Error! File '" + filename + "' is currently open for write. Cannot share with other output managers.\n");
      register_file(filename,Write);

      // Note: time has an unknown length. Setting its "length" to 0 tells the scorpio to
      // set this dimension as having an 'unlimited' length, thus allowing us to write
      // as many timesnaps to file as we desire.
      register_dimension(filename,"time","time",0);

      // Register time as a variable.
      register_variable(filename,"time","time",1,{"time"},  PIO_REAL,"time");

      // Make all output streams register their dims/vars
      for (auto& it : m_output_streams) {
        it->setup_output_file(filename);
      }

      // Set degree of freedom for "time"
      int time_dof[1] = {0};
      set_dof(filename,"time",0,time_dof);

      // Finish the definition phase for this file.
      eam_pio_enddef (filename); 
      if (is_checkpoint_step) { 
        set_int_attribute_c2f (filename.c_str(),"avg_count",m_output_control.nsteps_since_last_write);
      }
      filespecs.is_open = true;
    }

    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (m_is_model_restart_output || is_checkpoint_step) {
      if (m_io_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::out | std::ofstream::trunc);  // Open rpointer file and clear contents
        rpointer << filename << std::endl;
      }
    }

    // Update time in the output file
    pio_update_time(filename,timestamp.seconds_from(m_t0));
  }

  // Run the output streams
  for (auto& it : m_output_streams) {
    // Note: filename might referencing an invalid string, but it's only used
    //       in case is_write_step=true, in which case it will *for sure* contain
    //       a valid file name.
    it->run(filename,is_write_step,m_output_control.nsteps_since_last_write);
  }

  if (is_write_step) {
    // We're adding one snapshot to the file
    ++filespecs.num_snapshots_in_file;

    // If we wrote an output, the checkpoint counter needs to be reset
    if (is_output_step) {
      m_checkpoint_control.nsteps_since_last_write = 0;
    }

    // Finish up any updates to output file
    sync_outfile(filename);

    // Check if we need to close the output file
    if (filespecs.file_is_full()) {
      eam_pio_closefile(filename);
      filespecs.num_snapshots_in_file = 0;
      filespecs.is_open = false;
      control.nsteps_since_last_write = 0;
    }
  }
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  // 
  // Swapping with an empty mgr is the easiest way to cleanup.
  OutputManager other;
  std::swap(*this,other);
}

std::string OutputManager::
compute_filename_root (const IOControl& control) const
{
  return m_casename + "." +
         m_avg_type + "." +
         control.frequency_units+ "_x" +
         std::to_string(control.frequency);
}

void OutputManager::
set_params (const ekat::ParameterList& params,
            const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs)
{
  m_params = params;
  m_is_model_restart_output = m_params.get<bool>("Model Restart",false);

  if (m_is_model_restart_output) {
    using vos_t = std::vector<std::string>;

    // We build some restart parameters internally
    m_avg_type = ekat::upper_case(m_params.get<std::string>("Averaging Type","Instant"));
    EKAT_REQUIRE_MSG (m_avg_type=="INSTANT",
        "Error! For restart output, the averaging type must be 'Instant'.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");
    m_output_file_specs.max_snapshots_in_file = m_params.get("Max Snapshots Per File",1);
    EKAT_REQUIRE_MSG (m_output_file_specs.max_snapshots_in_file==1,
        "Error! For restart output, max snapshots per file must be 1.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");

    auto& fields_pl = m_params.sublist("Fields");
    for (const auto& it : field_mgrs) {
      auto restart_group = it.second->get_groups_info().at("RESTART");
      auto& fnames = fields_pl.get<vos_t>(it.first);
      EKAT_REQUIRE_MSG (fnames.size()==0,
        "Error! For restart output, don't specify the fields names. We can create this info internally.\n");
      for (const auto& n : restart_group->m_fields_names) {
        fnames.push_back(n);
      }
    }
    m_casename = m_params.get<std::string>("Casename", "scream_restart");
  } else {
    m_avg_type = ekat::upper_case(m_params.get<std::string>("Averaging Type"));
    m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File");
    m_casename = m_params.get<std::string>("Casename");
  }
}
/*===============================================================================================*/
} // namespace scream
