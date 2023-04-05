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
  m_output_control.frequency  = out_control_pl.get<int>("Frequency");
  m_output_control.frequency_units = out_control_pl.get<std::string>("frequency_units");
  m_output_control.timestamp_of_last_write   = start_ref ? m_case_t0 : m_run_t0;

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
    std::set<std::shared_ptr<const AbstractGrid>> grids;
    for (auto& it : m_output_streams) {
      grids.insert(it->get_io_grid());
    }

    // If 2+ grids are present, we mandate suffix on all geo_data fields,
    // to avoid clashes of names.
    bool use_suffix = grids.size()>1;
    for (auto grid : grids) {
      std::vector<Field> fields;
      for (const auto& fn : grid->get_geometry_data_names()) {
        const auto& f = grid->get_geometry_data(fn);
        if (use_suffix) {
          fields.push_back(f.clone(f.name()+"_"+grid->m_short_name));
        } else {
          fields.push_back(f.clone());
        }
      }
      auto output = std::make_shared<output_type>(m_io_comm,fields,grid);
      m_geo_data_streams.push_back(output);
    }
  }

  const auto has_restart_data = (m_avg_type!=OutputAvgType::Instant && m_output_control.frequency>1);
  if (has_restart_data && m_params.isSublist("Checkpoint Control")) {
    // Output control
    // TODO: It would be great if there was an option where, if Checkpoint Control was not a sublist, we
    //       could query the restart control information and just use that. 
    auto& pl = m_params.sublist("Checkpoint Control");
    m_checkpoint_control.frequency                 = pl.get<int>("Frequency");
    m_checkpoint_control.frequency_units           = pl.get<std::string>("frequency_units");
    m_checkpoint_control.timestamp_of_last_write   = run_t0;

    // File specs
    m_checkpoint_file_specs.max_snapshots_in_file = 1;
    m_checkpoint_file_specs.filename_with_mpiranks    = pl.get("MPI Ranks in Filename",false);
    m_checkpoint_file_specs.save_grid_data = false;
  } else {
    // If there is no restart data or there is but no checkpoint control sublist then we initialize
    // the checkpoint control so that it never writes checkpoints.
    m_checkpoint_control.frequency_units = "never";
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

    if (perform_history_restart) {
      // We can use the step counter in run_t0 to check at what point within an output interval
      // the previous simulation was stopped at.
      // NOTE: if you change the output frequency when you restart, this could lead to wonky behavior
      m_output_control.nsamples_since_last_write = m_run_t0.get_num_steps() % m_output_control.frequency; 

      if (has_restart_data) {
        using namespace scorpio;
        auto fn = find_filename_in_rpointer(hist_restart_casename,false,m_io_comm,m_run_t0);
        register_file(fn.c_str(),FileMode::Read);
        auto date = get_int_attribute_c2f (fn.c_str(), "last_write_date");
        auto time = get_int_attribute_c2f (fn.c_str(), "last_write_time");

        std::vector<int> vdate = {date/10000, (date/100)%100, date%100};
        std::vector<int> vtime = {time/10000, (time/100)%100, time%100};
        util::TimeStamp last_write_ts (vdate,vtime);
        m_output_control.timestamp_of_last_write = last_write_ts;
        eam_pio_closefile(fn.c_str());

        // If the type/freq of output needs restart data, we need to restart the streams
        if (m_output_control.nsamples_since_last_write>0) {
          for (auto stream : m_output_streams) {
            stream->restart(fn);
          }
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

void OutputManager::setup_globals_map (const globals_map_t& globals) {
  EKAT_REQUIRE_MSG (m_globals.size()==0,
      "Error! Globals already set in this output manager.\n");

  m_globals = globals;
}

/*===============================================================================================*/
void OutputManager::run(const util::TimeStamp& timestamp)
{
  using namespace scorpio;

  std::string timer_root = m_is_model_restart_output ? "EAMxx::IO::restart" : "EAMxx::IO::standard";
  start_timer(timer_root); 
  // Check if we need to open a new file
  ++m_output_control.nsamples_since_last_write;
  ++m_checkpoint_control.nsamples_since_last_write;

  // Check if this is a write step (and what kind)
  const bool is_t0_output       = timestamp==m_case_t0;
  const bool is_output_step     = m_output_control.is_write_step(timestamp) || is_t0_output;
  const bool is_checkpoint_step = m_checkpoint_control.is_write_step(timestamp) && not is_output_step;
  const bool is_write_step      = is_output_step || is_checkpoint_step;

  // If neither output or checkpoint, these won't be used anyways,
  // so no need to check if is_write_step == true.
  auto& control   = is_checkpoint_step ? m_checkpoint_control : m_output_control;
  auto& filespecs = is_checkpoint_step ? m_checkpoint_file_specs : m_output_file_specs;
  auto& filename  = filespecs.filename;

  // Compute filename (if write step)
  start_timer(timer_root+"::get_new_file"); 
  if (is_write_step) {
    // Check if we need to open a new file
    if (not filespecs.is_open) {
      // Register all dims/vars, write geometry data (e.g. lat/lon/hyam/hybm)
      setup_file(filespecs,control,timestamp);
    }

    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (m_is_model_restart_output || is_checkpoint_step) {
      if (m_io_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append to it
        rpointer << filename << std::endl;
      }
    }

    // Update time and nsteps in the output file
    pio_update_time(filename,timestamp.days_from(m_case_t0));
    if (m_is_model_restart_output) {
      // Only write nsteps on model restart
      set_int_attribute_c2f(filename.c_str(),"nsteps",timestamp.get_num_steps());
    } else if (is_checkpoint_step) {
      // Update the date of last write
      const auto& last_write_ts = control.timestamp_of_last_write;
      auto last_write_date = last_write_ts.get_date()[0]*10000 + last_write_ts.get_date()[1]*100 + last_write_ts.get_date()[2];
      auto last_write_time = last_write_ts.get_time()[0]*10000 + last_write_ts.get_time()[1]*100 + last_write_ts.get_time()[2];
      set_int_attribute_c2f(filename.c_str(),"last_write_date",last_write_date);
      set_int_attribute_c2f(filename.c_str(),"last_write_time",last_write_time);
    }
  }
  stop_timer(timer_root+"::get_new_file"); 

  // Log if we write output this step:
  if (m_atm_logger && is_write_step) {
    m_atm_logger->info("[EAMxx::output_manager] - Writing output:");
    m_atm_logger->info("[EAMxx::output_manager]      CASE: " + m_casename); 
    m_atm_logger->info("[EAMxx::output_manager]      FILE: " + filename); 
  }

  // Run the output streams
  start_timer(timer_root+"::run_output_streams"); 
  for (auto& it : m_output_streams) {
    // Note: filename might reference an invalid string, but it's only used
    //       in case is_write_step=true, in which case it will *for sure* contain
    //       a valid file name.
    it->run(filename,is_write_step,m_output_control.nsamples_since_last_write,is_t0_output);
  }
  stop_timer(timer_root+"::run_output_streams"); 

  if (is_write_step) {
    for (const auto& it : m_globals) {
      const auto& name = it.first;
      const auto& type_any = it.second;
      const auto& type = type_any.first;
      const auto& any = type_any.second;
      if (type=="int") {
        const int& value = ekat::any_cast<int>(any);
        set_int_attribute_c2f(filename.c_str(),name.c_str(),value);
      } else {
        EKAT_ERROR_MSG ("Error! Unsupported global attribute type.\n"
            " - file name  : " + filename + "\n"
            " - global name: " + name + "'\n"
            " - global type: " + type + "'\n");
      }
    }

    start_timer(timer_root+"::update_snapshot_tally"); 
    // We're adding one snapshot to the file
    ++filespecs.num_snapshots_in_file;

    if (m_time_bnds.size()>0) {
      m_time_bnds[1] = timestamp.days_from(m_case_t0);
      scorpio::grid_write_data_array(filename, "time_bnds", m_time_bnds.data(), 2);
      m_time_bnds[0] = m_time_bnds[1];
    }

    // Since we wrote to file we need to reset the nsamples_since_last_write, the timestamp ...
    control.nsamples_since_last_write = 0;
    control.timestamp_of_last_write = timestamp;
    // ... and the local views - unless it is a checkpoint step, then we keep local views.
    if (!is_checkpoint_step) {
      for (auto& it : m_output_streams) {
        it->reset_dev_views();
      }
    }

    // Check if we need to close the output file
    if (filespecs.file_is_full()) {
      eam_pio_closefile(filename);
      filespecs.num_snapshots_in_file = 0;
      filespecs.is_open = false;
    }

    // Whether we wrote an output or a checkpoint, the checkpoint counter needs to be reset
    m_checkpoint_control.nsamples_since_last_write = 0;
    m_checkpoint_control.timestamp_of_last_write = timestamp;

    stop_timer(timer_root+"::update_snapshot_tally"); 
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
  // NOTE: we do *not* allow this for checkpoints, since it would be risky if it gets somehow enabled
  //       inside an ERP cime test.
  if (not is_checkpoint_step && file_specs.filename_with_mpiranks) {
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
  set_int_attribute_c2f(filename.c_str(),"start_date",t0_date);
  set_int_attribute_c2f(filename.c_str(),"start_time",t0_time);
  set_str_attribute_c2f(filename.c_str(),"averaging_type",e2str(m_avg_type).c_str());
  set_str_attribute_c2f(filename.c_str(),"averaging_frequency_units",m_output_control.frequency_units.c_str());
  set_int_attribute_c2f(filename.c_str(),"averaging_frequency",m_output_control.frequency);
  set_int_attribute_c2f(filename.c_str(),"max_snapshots_per_file",m_output_file_specs.max_snapshots_in_file);
  set_file_header(filename);
  eam_pio_enddef (filename); 

  if (m_avg_type!=OutputAvgType::Instant) {
    // Unfortunately, attributes cannot be set in define mode (why?), so this could
    // not be done while we were setting the time_bnds
    set_int_attribute_c2f(filename.c_str(),"sample_size",control.frequency);
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
  set_str_attribute_c2f(filename.c_str(),"source","E3SM Atmosphere Model Version 4 (EAMxx)");  // TODO: probably want to make sure that new versions are reflected here.
  set_str_attribute_c2f(filename.c_str(),"case","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"title","EAMxx History File");
  set_str_attribute_c2f(filename.c_str(),"compset","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"git_hash","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"host","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"version","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"initial_file","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"topography_file","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"contact","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"institution_id","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"product","");  // TODO
  set_str_attribute_c2f(filename.c_str(),"component","ATM");
  set_str_attribute_c2f(filename.c_str(),"conventions","");  // TODO

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
