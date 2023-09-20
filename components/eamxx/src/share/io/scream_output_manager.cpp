#include "scream_output_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <fstream>
#include <memory>

namespace scream
{

// Local helper functions:
void set_file_header(const std::string& filename);

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& run_t0,
       const util::TimeStamp& case_t0,
       const bool is_model_restart_output)
{
  using map_t = std::map<std::string,std::shared_ptr<fm_type>>;
  map_t fms;
  fms[field_mgr->get_grid()->name()] = field_mgr;
  setup(io_comm,params,fms,grids_mgr,run_t0,case_t0,is_model_restart_output);
}

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& run_t0,
       const util::TimeStamp& case_t0,
       const bool is_model_restart_output)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (run_t0.is_valid(),
      "Error! Invalid case_t0 timestamp: " + case_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (run_t0.is_valid(),
      "Error! Invalid run_t0 timestamp: " + case_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (case_t0<=run_t0,
      "Error! The case_t0 timestamp must precede run_t0.\n"
      "   run_t0 : " + run_t0.to_string() + "\n"
      "   case_t0: " + case_t0.to_string() + "\n");

  m_io_comm = io_comm;
  m_run_t0 = run_t0;
  m_case_t0 = case_t0;
  m_is_restarted_run = (case_t0<run_t0);
  m_is_model_restart_output = is_model_restart_output;

  // Check for model restart output
  set_params(params,field_mgrs);

  // Output control
  EKAT_REQUIRE_MSG(m_params.isSublist("output_control"),
      "Error! The output control YAML file for " + m_casename + " is missing the sublist 'output_control'");
  auto& out_control_pl = m_params.sublist("output_control");
  // Determine which timestamp to use a reference for output frequency.  Two options:
  // 	1. use_case_as_start_reference: TRUE  - implies we want to calculate frequency from the beginning of the whole simulation, even if this is a restarted run.
  // 	2. use_case_as_start_reference: FALSE - implies we want to base the frequency of output on when this particular simulation started.
  // Note, (2) is needed for restarts since the restart frequency in CIME assumes a reference of when this run began.
  const bool start_ref = out_control_pl.get<bool>("use_case_as_start_reference",!m_is_model_restart_output);
  m_output_control.frequency_units = out_control_pl.get<std::string>("frequency_units");
  // In case output is disabled, no point in doing anything else
  if (m_output_control.frequency_units=="none" || m_output_control.frequency_units=="never") {
    m_output_disabled = true;
    return;
  }
  m_output_control.frequency = out_control_pl.get<int>("Frequency");
  EKAT_REQUIRE_MSG (m_output_control.frequency>0,
      "Error! Invalid frequency (" + std::to_string(m_output_control.frequency) + ") in Output Control. Please, use positive number.\n");

  m_output_control.timestamp_of_last_write = start_ref ? m_case_t0 : m_run_t0;

  // File specs
  m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File",-1);
  m_output_file_specs.filename_with_mpiranks    = out_control_pl.get("MPI Ranks in Filename",false);
  m_output_file_specs.save_grid_data            = out_control_pl.get("save_grid_data",!m_is_model_restart_output);

  // For each grid, create a separate output stream.
  if (field_mgrs.size()==1) {
    auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.begin()->second,grids_mgr);
    m_output_streams.push_back(output);
  } else {
    const auto& fields_pl = m_params.sublist("Fields");
    for (auto it=fields_pl.sublists_names_cbegin(); it!=fields_pl.sublists_names_cend(); ++it) {
      const auto& gname = *it;
      EKAT_REQUIRE_MSG (grids_mgr->has_grid(gname),
          "Error! Output requested on grid '" + gname + "', but the grids manager does not store such grid.\n");

      EKAT_REQUIRE_MSG (field_mgrs.find(gname)!=field_mgrs.end(),
          "Error! Output requested on grid '" + gname + "', but no field manager is available for such grid.\n");

      auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.at(gname),grids_mgr);
      m_output_streams.push_back(output);
    }
  }

  // For normal output, setup the geometry data streams, which we used to write the
  // geo data in the output file when we create it.
  if (m_output_file_specs.save_grid_data) {
    std::map<std::string, std::shared_ptr<const AbstractGrid>> grids;
    for (const auto& it : m_output_streams) {
      grids[it->get_io_grid()->name()] = it->get_io_grid();
    }

    // If 2+ grids are present, we mandate suffix on all geo_data fields,
    // to avoid clashes of names.
    bool use_suffix = grids.size()>1;
    for (const auto& grid : grids) {
      std::vector<Field> fields;
      for (const auto& fn : grid.second->get_geometry_data_names()) {
        const auto& f = grid.second->get_geometry_data(fn);
        if (use_suffix) {
          fields.push_back(f.clone(f.name()+"_"+grid.second->m_short_name));
        } else {
          fields.push_back(f.clone());
        }
      }
      auto output = std::make_shared<output_type>(m_io_comm,fields,grid.second);
      m_geo_data_streams.push_back(output);
    }
  }

  if (m_params.isSublist("Checkpoint Control")) {
    // Output control
    // TODO: It would be great if there was an option where, if Checkpoint Control was not a sublist, we
    //       could query the restart control information and just use that.
    auto& pl = m_params.sublist("Checkpoint Control");
    m_checkpoint_control.frequency_units           = pl.get<std::string>("frequency_units");

    if (m_checkpoint_control.output_enabled()) {
      m_checkpoint_control.timestamp_of_last_write   = run_t0;
      m_checkpoint_control.frequency = pl.get<int>("Frequency");
      EKAT_REQUIRE_MSG (m_output_control.frequency>0,
          "Error! Invalid frequency (" + std::to_string(m_checkpoint_control.frequency) + ") in Checkpoint Control. Please, use positive number.\n");

      // File specs
      m_checkpoint_file_specs.max_snapshots_in_file = 1;
      m_checkpoint_file_specs.filename_with_mpiranks    = pl.get("MPI Ranks in Filename",false);
      m_checkpoint_file_specs.save_grid_data = false;
      m_checkpoint_file_specs.hist_restart_file = true;
    }
  }

  // If this is normal output (not the model restart output) and the output specs
  // require it, we need to restart the output history.
  // E.g., we might save 30-day avg value for field F, but due to job size
  // break the run into three 10-day runs. We then need to save the state of
  // our averaging in a "restart" file (e.g., the current avg).
  // Note: the user might decide *not* to restart the output, so give the option
  //       of disabling the restart. Also, the user might want to change the
  //       filename_prefix, so allow to specify a different filename_prefix for the restart file.
  if (m_is_restarted_run) {
    // Allow to skip history restart, or to specify a filename_prefix for the restart file
    // that is different from the filename_prefix of the current output.
    auto& restart_pl = m_params.sublist("Restart");
    bool perform_history_restart = restart_pl.get("Perform Restart",true);
    auto hist_restart_casename = restart_pl.get("filename_prefix",m_casename);

    if (m_is_model_restart_output) {
      // For model restart output, the restart time (which is the start time of this run) is precisely
      // when the last write happened, so we can quickly init the output control.
      m_output_control.timestamp_of_last_write = m_run_t0;
      m_output_control.nsamples_since_last_write = 0;
    } else if (perform_history_restart) {
      using namespace scorpio;
      auto fn = find_filename_in_rpointer(hist_restart_casename,false,m_io_comm,m_run_t0);

      // From restart file, get the time of last write, as well as the current size of the avg sample
      m_output_control.timestamp_of_last_write = read_timestamp(fn,"last_write");
      m_output_control.nsamples_since_last_write = get_attribute<int>(fn,"num_snapshots_since_last_write");

      // If the type/freq of output needs restart data, we need to restart the streams
      const auto has_restart_data = m_avg_type!=OutputAvgType::Instant && m_output_control.frequency>1;
      if (has_restart_data && m_output_control.nsamples_since_last_write>0) {
        for (auto stream : m_output_streams) {
          stream->restart(fn);
        }
      }
    }
  }

  if (m_avg_type!=OutputAvgType::Instant) {
    // Init the left hand point of time_bnds based on run/case t0.
    m_time_bnds.resize(2);
    m_time_bnds[0] = m_run_t0.days_from(m_case_t0);
  } else if (m_output_control.output_enabled() && m_run_t0==m_case_t0 && !m_is_model_restart_output) {
    this->run(m_run_t0);
  }

  // Log this output stream
  push_to_logger();
}

void OutputManager::
add_global (const std::string& name, const ekat::any& global) {
  EKAT_REQUIRE_MSG (m_globals.find(name)==m_globals.end(),
      "Error! Global attribute was already set in this output manager.\n"
      "  - global att name: " + name + "\n");

  m_globals[name] = global;
}

/*===============================================================================================*/
void OutputManager::run(const util::TimeStamp& timestamp)
{
  // In case output is disabled, no point in doing anything else
  if (m_output_disabled) {
    return;
  }
  using namespace scorpio;

  std::string timer_root = m_is_model_restart_output ? "EAMxx::IO::restart" : "EAMxx::IO::standard";
  start_timer(timer_root);
  // Check if we need to open a new file
  ++m_output_control.nsamples_since_last_write;
  ++m_checkpoint_control.nsamples_since_last_write;

  // Check if this is a write step (and what kind)
  // Note: a full checkpoint not only writes globals in the restart file, but also all the history variables.
  //       Since we *always* write a history restart file, we can have a non-full checkpoint, if the average
  //       type is Instant and/or the frequency is every step. A non-full checkpoint will simply write some
  //       global attribute, such as the time of last write.
  const bool is_t0_output            = timestamp==m_case_t0;
  const bool is_output_step          = m_output_control.is_write_step(timestamp) || is_t0_output;
  const bool is_checkpoint_step      = m_checkpoint_control.is_write_step(timestamp);
  const bool has_checkpoint_data     = (m_avg_type!=OutputAvgType::Instant && m_output_control.frequency>1);
  const bool is_full_checkpoint_step = is_checkpoint_step && has_checkpoint_data && not is_output_step;
  const bool is_write_step           = is_output_step || is_checkpoint_step;

  // Create and setup output/checkpoint file(s), if necessary
  start_timer(timer_root+"::get_new_file");
  auto setup_output_file = [&](IOControl& control, IOFileSpecs& filespecs, bool add_to_rpointer, const std::string& file_type) {
    // Check if we need to open a new file
    if (not filespecs.is_open) {
      // Register all dims/vars, write geometry data (e.g. lat/lon/hyam/hybm)
      setup_file(filespecs,control,timestamp);
    }

    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (add_to_rpointer) {
      if (m_io_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append to it
        rpointer << filespecs.filename << std::endl;
      }
    }

    if (m_atm_logger) {
      m_atm_logger->info("[EAMxx::output_manager] - Writing " + file_type + ":");
      m_atm_logger->info("[EAMxx::output_manager]      CASE: " + m_casename);
      m_atm_logger->info("[EAMxx::output_manager]      FILE: " + filespecs.filename);
    }
  };

  if (is_output_step) {
    setup_output_file(m_output_control,m_output_file_specs,m_is_model_restart_output,m_is_model_restart_output ? "model restart" : "model output");

    // Update time (must be done _before_ writing fields)
    pio_update_time(m_output_file_specs.filename,timestamp.days_from(m_case_t0));
  }
  if (is_checkpoint_step) {
    setup_output_file(m_checkpoint_control,m_checkpoint_file_specs,true,"history restart");

    if (is_full_checkpoint_step) {
      // Update time (must be done _before_ writing fields)
      pio_update_time(m_checkpoint_file_specs.filename,timestamp.days_from(m_case_t0));
    }
  }
  stop_timer(timer_root+"::get_new_file");

  // Run the output streams
  start_timer(timer_root+"::run_output_streams");
  const auto& fields_write_filename = is_output_step ? m_output_file_specs.filename : m_checkpoint_file_specs.filename;
  for (auto& it : m_output_streams) {
    // Note: filename only matters if is_output_step || is_full_checkpoint_step=true. In that case, it will definitely point to a valid file name.
    it->run(fields_write_filename,is_output_step || is_full_checkpoint_step,m_output_control.nsamples_since_last_write,is_t0_output);
  }
  stop_timer(timer_root+"::run_output_streams");

  if (is_write_step) {
    if (m_time_bnds.size()>0) {
      m_time_bnds[1] = timestamp.days_from(m_case_t0);
    }

    // If we write output, reset local views
    if (is_output_step) {
      for (auto& it : m_output_streams) {
        it->reset_dev_views();
      }
    }

    auto write_global_data = [&](IOControl& control, IOFileSpecs& filespecs) {
      if (m_is_model_restart_output) {
        // Only write nsteps on model restart
        set_attribute(filespecs.filename,"nsteps",timestamp.get_num_steps());
      } else if (filespecs.hist_restart_file) {
        // Update the date of last write and sample size
        scorpio::write_timestamp (filespecs.filename,"last_write",m_output_control.timestamp_of_last_write);
        scorpio::set_attribute (filespecs.filename,"num_snapshots_since_last_write",m_output_control.nsamples_since_last_write);
      }

      // Write all stored globals
      for (const auto& it : m_globals) {
        const auto& name = it.first;
        const auto& any = it.second;
        set_any_attribute(filespecs.filename,name,any);
      }

      // We're adding one snapshot to the file
      ++filespecs.num_snapshots_in_file;

      if (m_time_bnds.size()>0) {
        scorpio::grid_write_data_array(filespecs.filename, "time_bnds", m_time_bnds.data(), 2);
      }

      // Since we wrote to file we need to reset the nsamples_since_last_write, the timestamp ...
      control.nsamples_since_last_write = 0;
      control.timestamp_of_last_write = timestamp;

      // Check if we need to close the output file
      if (filespecs.file_is_full()) {
        eam_pio_closefile(filespecs.filename);
        filespecs.num_snapshots_in_file = 0;
        filespecs.is_open = false;
      }
    };

    start_timer(timer_root+"::update_snapshot_tally");
    // Important! Process output file first, and hist restart (if any) second.
    // That's b/c write_global_data will update m_output_control.timestamp_of_last_write,
    // which is later be written as global data in the hist restart file
    if (is_output_step) {
      write_global_data(m_output_control,m_output_file_specs);
    }
    if (is_checkpoint_step) {
      write_global_data(m_checkpoint_control,m_checkpoint_file_specs);
    }
    stop_timer(timer_root+"::update_snapshot_tally");
    if (m_time_bnds.size()>0) {
      m_time_bnds[0] = m_time_bnds[1];
    }
  }

  stop_timer(timer_root);
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  // Close any output file still open
  if (m_output_file_specs.is_open) {
    scorpio::eam_pio_closefile (m_output_file_specs.filename);
  }
  if (m_checkpoint_file_specs.is_open) {
    scorpio::eam_pio_closefile (m_checkpoint_file_specs.filename);
  }

  // Swapping with an empty mgr is the easiest way to cleanup.
  OutputManager other;
  std::swap(*this,other);
}

long long OutputManager::res_dep_memory_footprint () const {
  long long mf = 0;
  for (const auto& os : m_output_streams) {
    mf += os->res_dep_memory_footprint();
  }

  return mf;
}

std::string OutputManager::
compute_filename (const IOControl& control,
                  const IOFileSpecs& file_specs,
                  const bool is_checkpoint_step,
                  const util::TimeStamp& timestamp) const
{
  std::string suffix =
    is_checkpoint_step ? ".rhist"
                       : (m_is_model_restart_output ? ".r" : "");
  auto filename = m_casename + suffix;

  // Always add avg type and frequency info
  filename += "." + e2str(m_avg_type);
  filename += "." + control.frequency_units+ "_x" + std::to_string(control.frequency);

  // Optionally, add number of mpi ranks (useful mostly in unit tests, to run multiple MPI configs in parallel)
  if (file_specs.filename_with_mpiranks) {
    filename += ".np" + std::to_string(m_io_comm.size());
  }

  // Always add a time stamp
  if (m_avg_type==OutputAvgType::Instant || is_checkpoint_step) {
    filename += "." + timestamp.to_string();
  } else {
    filename += "." + control.timestamp_of_last_write.to_string();
  }

  return filename + ".nc";
}

void OutputManager::
set_params (const ekat::ParameterList& params,
            const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs)
{
  m_params = params;
  if (m_is_model_restart_output) {
    using vos_t = std::vector<std::string>;

    // We build some restart parameters internally
    auto avg_type = m_params.get<std::string>("Averaging Type","INSTANT");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type==OutputAvgType::Instant,
        "Error! For restart output, the averaging type must be 'Instant'.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");

    m_output_file_specs.max_snapshots_in_file = m_params.get("Max Snapshots Per File",1);
    EKAT_REQUIRE_MSG (m_output_file_specs.max_snapshots_in_file==1,
        "Error! For restart output, max snapshots per file must be 1.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");

    auto& fields_pl = m_params.sublist("Fields");
    for (const auto& it : field_mgrs) {
      const auto& fm = it.second;
      vos_t fnames;
      // There may be no RESTART group on this grid
      if (fm->has_group("RESTART")) {
        auto restart_group = fm->get_groups_info().at("RESTART");
        EKAT_REQUIRE_MSG (not fields_pl.isParameter(it.first),
          "Error! For restart output, don't specify the fields names. We will create this info internally.\n");
        for (const auto& n : restart_group->m_fields_names) {
          fnames.push_back(n);
        }
      }
      fields_pl.sublist(it.first).set("Field Names",fnames);
    }
    m_casename = m_params.get<std::string>("filename_prefix");
    // Match precision of Fields
    m_params.set<std::string>("Floating Point Precision","real");
  } else {
    auto avg_type = m_params.get<std::string>("Averaging Type");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
        "Error! Unsupported averaging type '" + avg_type + "'.\n"
        "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

    m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File",-1);
    m_casename = m_params.get<std::string>("filename_prefix");

    // Allow user to ask for higher precision for normal model output,
    // but default to single to save on storage
    const auto& prec = m_params.get<std::string>("Floating Point Precision", "single");
    EKAT_REQUIRE_MSG (prec=="single" || prec=="double" || prec=="real",
        "Error! Invalid floating point precision '" + prec + "'.\n");
  }
}
/*===============================================================================================*/
void OutputManager::
setup_file (      IOFileSpecs& filespecs, const IOControl& control,
            const util::TimeStamp& timestamp)
{
  using namespace scorpio;

  const bool is_checkpoint_step = &control==&m_checkpoint_control;
  auto& filename = filespecs.filename;

  // Compute new file name
  filename = compute_filename (control,filespecs,is_checkpoint_step,timestamp);

  // Register new netCDF file for output. First, check no other output managers
  // are trying to write on the same file
  EKAT_REQUIRE_MSG (not is_file_open_c2f(filename.c_str(),Write),
      "Error! File '" + filename + "' is currently open for write. Cannot share with other output managers.\n");
  register_file(filename,Write);

  // Note: time has an unknown length. Setting its "length" to 0 tells the scorpio to
  // set this dimension as having an 'unlimited' length, thus allowing us to write
  // as many timesnaps to file as we desire.
  register_dimension(filename,"time","time",0,false);

  // Register time as a variable.
  auto time_units="days since " + m_case_t0.get_date_string() + " " + m_case_t0.get_time_string();
  register_variable(filename,"time","time",time_units,{"time"}, "double", "double","time");
#ifdef SCREAM_HAS_LEAP_YEAR
  set_variable_metadata (filename,"time","calendar","gregorian");
#else
  set_variable_metadata (filename,"time","calendar","noleap");
#endif
  if (m_avg_type!=OutputAvgType::Instant) {
    // First, ensure a 'dim2' dimension with len=2 is registered.
    register_dimension(filename,"dim2","dim2",2,false);

    // Register time_bnds var, with its dofs
    register_variable(filename,"time_bnds","time_bnds",time_units,{"dim2","time"},"double","double","time-dim2");
    scorpio::offset_t time_bnds_dofs[2] = {0,1};
    set_dof(filename,"time_bnds",2,time_bnds_dofs);

    // Make it clear how the time_bnds should be interpreted
    set_variable_metadata(filename,"time_bnds","note","right endpoint accummulation");

    // I'm not sure what's the point of this, but CF conventions seem to require it
    set_variable_metadata (filename,"time","bounds","time_bnds");
  }

  std::string fp_precision = is_checkpoint_step
                           ? "real"
                           : m_params.get<std::string>("Floating Point Precision");

  // Make all output streams register their dims/vars
  for (auto& it : m_output_streams) {
    it->setup_output_file(filename,fp_precision);
  }

  if (filespecs.save_grid_data) {
    // If not a restart file, also register geo data fields.
    for (auto& it : m_geo_data_streams) {
      it->setup_output_file(filename,fp_precision);
    }
  }

  // Set degree of freedom for "time"
  scorpio::offset_t time_dof[1] = {0};
  set_dof(filename,"time",0,time_dof);

  // Finish the definition phase for this file.
  auto t0_date = m_case_t0.get_date()[0]*10000 + m_case_t0.get_date()[1]*100 + m_case_t0.get_date()[2];
  auto t0_time = m_case_t0.get_time()[0]*10000 + m_case_t0.get_time()[1]*100 + m_case_t0.get_time()[2];

  set_attribute(filename,"start_date",t0_date);
  set_attribute(filename,"start_time",t0_time);
  set_attribute(filename,"averaging_type",e2str(m_avg_type));
  set_attribute(filename,"averaging_frequency_units",m_output_control.frequency_units);
  set_attribute(filename,"averaging_frequency",m_output_control.frequency);
  set_attribute(filename,"max_snapshots_per_file",m_output_file_specs.max_snapshots_in_file);
  set_file_header(filename);
  eam_pio_enddef (filename);

  if (m_avg_type!=OutputAvgType::Instant) {
    // Unfortunately, attributes cannot be set in define mode (why?), so this could
    // not be done while we were setting the time_bnds
    set_attribute(filename,"sample_size",control.frequency);
  }

  if (filespecs.save_grid_data) {
    // Immediately run the geo data streams
    for (const auto& it : m_geo_data_streams) {
      it->run(filename,true,0);
    }
  }

  filespecs.is_open = true;
}
/*===============================================================================================*/
void set_file_header(const std::string& filename)
{
  using namespace scorpio;

  // TODO: All attributes marked TODO below need to be set.  Hopefully by a universal value that reflects
  // what the attribute is.  For example, git-hash should be the git-hash associated with this version of
  // the code at build time for this executable.
  set_attribute<std::string>(filename,"source","E3SM Atmosphere Model Version 4 (EAMxx)");  // TODO: probably want to make sure that new versions are reflected here.
  set_attribute<std::string>(filename,"case","");  // TODO
  set_attribute<std::string>(filename,"title","EAMxx History File");
  set_attribute<std::string>(filename,"compset","");  // TODO
  set_attribute<std::string>(filename,"git_hash","");  // TODO
  set_attribute<std::string>(filename,"host","");  // TODO
  set_attribute<std::string>(filename,"version","");  // TODO
  set_attribute<std::string>(filename,"initial_file","");  // TODO
  set_attribute<std::string>(filename,"topography_file","");  // TODO
  set_attribute<std::string>(filename,"contact","");  // TODO
  set_attribute<std::string>(filename,"institution_id","");  // TODO
  set_attribute<std::string>(filename,"product","");  // TODO
  set_attribute<std::string>(filename,"component","ATM");
  set_attribute<std::string>(filename,"Conventions","CF-1.8");  // TODO: In the future we may be able to have this be set at runtime.  We hard-code for now, because post-processing needs something in this global attribute. 2023-04-12
}
/*===============================================================================================*/
void OutputManager::
push_to_logger()
{
  // If no atm logger set then don't do anything
  if (!m_atm_logger) return;

  auto bool_to_string = [](const bool x) {
    std::string y = x ? "YES" : "NO";
    return y;
  };

  m_atm_logger->info("[EAMxx::output_manager] - New Output stream");
  m_atm_logger->info("                      Case: " + m_casename);
  m_atm_logger->info("                    Run t0: " + m_run_t0.to_string());
  m_atm_logger->info("                   Case t0: " + m_case_t0.to_string());
  m_atm_logger->info("              Reference t0: " + m_output_control.timestamp_of_last_write.to_string());
  m_atm_logger->info("         Is Restart File ?: " + bool_to_string(m_is_model_restart_output));
  m_atm_logger->info("        Is Restarted Run ?: " + bool_to_string(m_is_restarted_run));
  m_atm_logger->info("            Averaging Type: " + e2str(m_avg_type));
  m_atm_logger->info("          Output Frequency: " + std::to_string(m_output_control.frequency) + " " + m_output_control.frequency_units);
  m_atm_logger->info("         Max snaps in file: " + std::to_string(m_output_file_specs.max_snapshots_in_file));  // TODO: add "not set" if the value is -1
  m_atm_logger->info("      Includes Grid Data ?: " + bool_to_string(m_output_file_specs.save_grid_data));
  // List each GRID - TODO
  // List all FIELDS - TODO



}

} // namespace scream
