#include "scream_output_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"
#include "share/scream_config.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <fstream>
#include <memory>
#include <chrono>
#include <ctime>

namespace scream
{

OutputManager::
~OutputManager ()
{
  finalize();
}

void OutputManager::
initialize(const ekat::Comm& io_comm, const ekat::ParameterList& params,
           const util::TimeStamp& run_t0, const util::TimeStamp& case_t0,
           const bool is_model_restart_output, const RunType run_type)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (run_t0.is_valid(),
      "Error! Invalid run_t0 timestamp: " + run_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (case_t0.is_valid(),
      "Error! Invalid case_t0 timestamp: " + case_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (case_t0<=run_t0,
      "Error! The case_t0 timestamp must precede run_t0.\n"
      "   run_t0 : " + run_t0.to_string() + "\n"
      "   case_t0: " + case_t0.to_string() + "\n");

  m_io_comm = io_comm;
  m_params = params;
  m_run_t0 = run_t0;
  m_case_t0 = case_t0;
  m_run_type = run_type;
  m_is_model_restart_output = is_model_restart_output;
}

void OutputManager::
setup (const std::shared_ptr<fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr)
{
  using map_t = std::map<std::string,std::shared_ptr<fm_type>>;
  map_t fms;
  fms[field_mgr->get_grid()->name()] = field_mgr;
  setup(fms,grids_mgr);
}

void OutputManager::
setup (const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
       const std::shared_ptr<const gm_type>& grids_mgr)
{
  // Read input parameters and setup internal data
  setup_internals(field_mgrs);

  // Here, store if PG2 fields will be present in output streams.
  // Will be useful if multiple grids are defined (see below).
  bool pg2_grid_in_io_streams = false;
  const auto& fields_pl = m_params.sublist("Fields");
  for (auto it=fields_pl.sublists_names_cbegin(); it!=fields_pl.sublists_names_cend(); ++it) {
    if (*it == "Physics PG2") pg2_grid_in_io_streams = true;
  }

  // For each grid, create a separate output stream.
  if (field_mgrs.size()==1) {
    auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.begin()->second,grids_mgr);
    output->set_logger(m_atm_logger);
    m_output_streams.push_back(output);
  } else {
    for (auto it=fields_pl.sublists_names_cbegin(); it!=fields_pl.sublists_names_cend(); ++it) {
      const auto& gname = *it;

      // If this is a GLL grid (or IO Grid is GLL) and PG2 fields
      // were found above, we must reset the grid COL tag name to
      // be "ncol_d" to avoid conflicting lengths with ncol on
      // the PG2 grid.
      if (pg2_grid_in_io_streams) {
        const auto& grid_pl = fields_pl.sublist(gname);
        bool reset_ncol_naming = false;
        if (gname == "Physics GLL") reset_ncol_naming = true;
        if (grid_pl.isParameter("IO Grid Name")) {
          if (grid_pl.get<std::string>("IO Grid Name") == "Physics GLL") {
            reset_ncol_naming = true;
          }
        }
        if (reset_ncol_naming) {
          grids_mgr->
            get_grid_nonconst(grid_pl.get<std::string>("IO Grid Name"))->
              reset_field_tag_name(ShortFieldTagsNames::COL,"ncol_d");
	      }
      }

      EKAT_REQUIRE_MSG (grids_mgr->has_grid(gname),
          "Error! Output requested on grid '" + gname + "', but the grids manager does not store such grid.\n");

      EKAT_REQUIRE_MSG (field_mgrs.find(gname)!=field_mgrs.end(),
          "Error! Output requested on grid '" + gname + "', but no field manager is available for such grid.\n");

      auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.at(gname),grids_mgr);
      output->set_logger(m_atm_logger);
      m_output_streams.push_back(output);
    }
  }

  // For normal output, setup the geometry data streams, which we used to write the
  // geo data in the output file when we create it.
  if (m_save_grid_data) {
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

        if (f.rank()==0) {
          // Right now, this only happens for `dx_short`, a single scalar
          // coming from iop. Since that scalar is easily recomputed
          // upon restart, we can skip this
          // NOTE: without this, the code crash while attempting to get a
          //       pio decomp for this var, since there are no dimensions.
          //       Perhaps you can explore setting the var as a global att
          continue;
        }
        if (use_suffix) {
          fields.push_back(f.clone(f.name()+"_"+grid.second->m_short_name));
        } else {
          fields.push_back(f.clone());
        }
      }

      // See comment above for ncol naming with 2+ grids
      auto grid_nonconst = grid.second->clone(grid.first,true);
      if (grid.first == "Physics GLL" && pg2_grid_in_io_streams) {
        grid_nonconst->reset_field_tag_name(ShortFieldTagsNames::COL,"ncol_d");
      }

      auto output = std::make_shared<output_type>(m_io_comm,fields,grid_nonconst);
      m_geo_data_streams.push_back(output);
    }
  }

  // If this is model output (not model restart) we need to restart the output history.
  // For instant output, this just entails restarting timestamps and counters, while
  // for average output we also need to restore the accumulation state of the output fields
  // Note: the user might decide *not* to restart the output, so give the option
  //       of disabling the restart. Also, the user might want to change the
  //       filename_prefix, so allow to specify a different filename_prefix for the restart file.
  if (m_run_type==RunType::Restart and not m_is_model_restart_output) {
    // Allow to skip history restart, or to specify a filename_prefix for the restart file
    // that is different from the filename_prefix of the current output.
    auto& restart_pl = m_params.sublist("Restart");
    bool perform_history_restart = restart_pl.get("Perform Restart",true);
    auto hist_restart_filename_prefix = restart_pl.get("filename_prefix",m_filename_prefix);

    if (perform_history_restart) {
      using namespace scorpio;
      IOFileSpecs hist_restart_specs;
      hist_restart_specs.ftype = FileType::HistoryRestart;
      auto rhist_file = find_filename_in_rpointer(hist_restart_filename_prefix,false,m_io_comm,m_run_t0,m_avg_type,m_output_control);

      scorpio::register_file(rhist_file,scorpio::Read);
      // From restart file, get the time of last write, as well as the current size of the avg sample
      m_output_control.last_write_ts = read_timestamp(rhist_file,"last_write",true);
      m_output_control.compute_next_write_ts();
      m_output_control.nsamples_since_last_write = get_attribute<int>(rhist_file,"GLOBAL","num_snapshots_since_last_write");

      if (m_avg_type!=OutputAvgType::Instant) {
        m_time_bnds.resize(2);
        m_time_bnds[0] = m_output_control.last_write_ts.days_from(m_case_t0);
      }

      // If the type/freq of output needs restart data, we need to restart the streams
      const bool output_every_step   = m_output_control.frequency_units=="nsteps" &&
                                       m_output_control.frequency==1;
      const bool has_checkpoint_data = m_avg_type!=OutputAvgType::Instant && not output_every_step;
      if (has_checkpoint_data && m_output_control.nsamples_since_last_write>0) {
        for (auto stream : m_output_streams) {
          stream->restart(rhist_file);
        }
      }

      // We do NOT allow changing output specs across restart. If you do want to change
      // any of these, you MUST start a new output stream (e.g., setting 'Perform Restart: false')
      // NOTE: we do not check that freq/freq_units/avg_type are not changed: since we used
      //       that info to find the correct rhist file, we already know that they match!
      auto old_storage_type = scorpio::get_attribute<std::string>(rhist_file,"GLOBAL","file_max_storage_type");
      EKAT_REQUIRE_MSG (old_storage_type == e2str(m_output_file_specs.storage.type),
          "Error! Cannot change file storage type when performing history restart.\n"
          "  - old file_max_storage_type: " << old_storage_type << "\n"
          "  - new file_max_storage_type: " << e2str(m_output_file_specs.storage.type) << "\n"
          "If you *really* want to change the file storage type, you need to force using a new file, setting\n"
          "  Restart:\n"
          "    force_new_file: true\n");
      if (old_storage_type=="num_snapshot") {
        auto old_max_snaps = scorpio::get_attribute<int>(rhist_file,"GLOBAL","max_snapshots_per_file");
        EKAT_REQUIRE_MSG (old_max_snaps == m_output_file_specs.storage.max_snapshots_in_file,
            "Error! Cannot change max snapshots per file when performing history restart.\n"
            "  - old max snaps: " << old_max_snaps << "\n"
            "  - new max snaps: " << m_output_file_specs.storage.max_snapshots_in_file << "\n"
            "If you *really* want to change the file capacity, you need to force using a new file, setting\n"
            "  Restart:\n"
            "    force_new_file: true\n");
      }
      std::string fp_precision = m_params.get<std::string>("Floating Point Precision");
      auto old_fp_precision = scorpio::get_attribute<std::string>(rhist_file,"GLOBAL","fp_precision");
      EKAT_REQUIRE_MSG (old_fp_precision == fp_precision,
          "Error! Cannot change floating point precision when performing history restart.\n"
          "  - old fp precision: " << old_fp_precision << "\n"
          "  - new fp precision: " << fp_precision << "\n");

      // Check if the prev run wrote any output file (it may have not, if the restart was written
      // before the 1st output step). If there is a file, check if there's still room in it.
      const auto& last_output_filename = get_attribute<std::string>(rhist_file,"GLOBAL","last_output_filename");
      m_resume_output_file = last_output_filename!="" and not restart_pl.get("force_new_file",false);
      if (m_resume_output_file) {
        m_output_file_specs.storage.num_snapshots_in_file = scorpio::get_attribute<int>(rhist_file,"GLOBAL","last_output_file_num_snaps");

        if (m_output_file_specs.storage.snapshot_fits(m_output_control.next_write_ts)) {
          // The setup_file call will not register any new variable (the file is in Append mode,
          // so all dims/vars must already be in the file). However, it will register decompositions,
          // since those are a property of the run, not of the file.
          m_output_file_specs.filename = last_output_filename;
          m_output_file_specs.is_open = true;
          setup_file(m_output_file_specs,m_output_control);
        } else {
          m_output_file_specs.close();
        }
      }
      scorpio::release_file(rhist_file);
    }
  }

  // If m_time_bnds.size()>0, it was already inited during restart
  if (m_avg_type!=OutputAvgType::Instant && m_time_bnds.size()==0) {
    // Init the left hand point of time_bnds based on run/case t0.
    m_time_bnds.resize(2);
    m_time_bnds[0] = m_run_t0.days_from(m_case_t0);
  } else if (m_output_control.output_enabled() and
             m_run_type==RunType::Initial and
             not m_is_model_restart_output and
             not m_params.sublist("output_control").get<bool>("skip_t0_output",false)) // This will be true for ERS/ERP tests
  {
    // In order to trigger a t0 write, we need to have next_write_ts matching run_t0
    m_output_control.next_write_ts = m_run_t0;
    // This is in case some diags need to init the timestep. Their output may be meaningless
    // at t0 (e.g., if their input fields are not in the initial condition fields set,
    // and have yet to be computed), but they may still require the start-of-step timestamp to be valid
    init_timestep(m_run_t0,0);
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
void OutputManager::init_timestep (const util::TimeStamp& start_of_step, const Real dt)
{
  // In case output is disabled, no point in doing anything else
  if (not m_output_control.output_enabled()) {
    return;
  }

  const bool is_first_step = m_output_control.dt==0 and dt>0;

  // Make sure dt is in the control
  m_output_control.set_dt(dt);

  if (m_run_type==RunType::Initial and is_first_step and m_avg_type==OutputAvgType::Instant and
      m_output_file_specs.storage.type!=NumSnaps and m_output_control.frequency_units=="nsteps") {
    // This is the 1st step of the whole run, and a very sneaky corner case. Bear with me.
    // When we call run, we also compute next_write_ts. Then, we use next_write_ts to see if the
    // next output step will fit in the currently open file, and, if not, close it right away.
    // For a storage type!=NumSnaps, we need to have a valid timestamp for next_write_ts, which
    // for freq=nsteps requires to know dt. But at t=case_t0, we did NOT have dt, which means we
    // computed next_write_ts=last_write_ts (in terms of date:time, the num_steps is correct).
    // This means that at that time we deemed that the next_write_ts definitely fit in the same
    // file as last_write_ts (date/time are the same!), which may or may not be true for non NumSnaps
    // storage. To fix this, we recompute next_write_ts here, and close the file if it doesn't.
    m_output_control.compute_next_write_ts();
    close_or_flush_if_needed (m_output_file_specs,m_output_control);
  }

  // Note: we need to "init" the timestep if we are going to do something this step, which means we either
  //       have INST output and it's a write step, or we have AVG output.
  const auto end_of_step = start_of_step+dt;
  const bool is_output_step = m_output_control.is_write_step(end_of_step) || end_of_step==m_case_t0;
  if (not is_output_step and m_avg_type==OutputAvgType::Instant) {
    return;
  }

  if (m_atm_logger) {
    m_atm_logger->debug("[OutputManager::init_timestep] filename_prefix: " + m_filename_prefix + "\n");
  }

  for (auto s : m_output_streams) {
    s->init_timestep(start_of_step);
  }
}

void OutputManager::run(const util::TimeStamp& timestamp)
{
  // In case output is disabled, no point in doing anything else
  if (not m_output_control.output_enabled()) {
    return;
  }

  // Ensure we did not go past the scheduled write time without hitting it
  EKAT_REQUIRE_MSG (
      (m_output_control.frequency_units=="nsteps"
          ? timestamp.get_num_steps()<=m_output_control.next_write_ts.get_num_steps()
          : timestamp<=m_output_control.next_write_ts),
      "Error! The input timestamp is past the next scheduled write timestamp.\n"
      "  - current time stamp   : " + timestamp.to_string() + "\n"
      "  - next write time stamp: " + m_output_control.next_write_ts.to_string() + "\n"
      "The most likely cause is an output frequency that is faster than the atm timestep.\n"
      "Try to increase 'Frequency' and/or 'frequency_units' in your output yaml file.\n");

  if (m_atm_logger) {
    m_atm_logger->debug("[OutputManager::run] filename_prefix: " + m_filename_prefix + "\n");
  }

  using namespace scorpio;

  std::string timer_root = m_is_model_restart_output ? "EAMxx::IO::restart" : "EAMxx::IO::standard";
  start_timer(timer_root);
  start_timer("EAMxx::IO::" + m_params.name());

  // Check if this is a write step (and what kind)
  // Note: a full checkpoint not only writes globals in the restart file, but also all the history variables.
  //       Since we *always* write a history restart file, we can have a non-full checkpoint, if the average
  //       type is Instant and/or the frequency is every step. A non-full checkpoint will simply write some
  //       global attribute, such as the time of last write.
  //       Also, notice that units="nhours" and freq=1 would still output evey step if dt=3600s. However,
  //       it is somewhat hard to figure out if output happens every step, without having a dt to compare
  //       against. Therefore, we simply assume that if units!=nsteps OR freq>1, then we don't output every
  //       timestep. If, in fact, we are outputing every timestep, it's likely a small test, so it's not too
  //       bad if we write out some extra data.
  const bool output_every_step       = m_output_control.frequency_units=="nsteps" &&
                                       m_output_control.frequency==1;
  const bool is_t0_output            = timestamp==m_case_t0;
  const bool is_output_step          = m_output_control.is_write_step(timestamp) || is_t0_output;
  const bool is_checkpoint_step      = m_checkpoint_control.is_write_step(timestamp) && not is_t0_output;
  const bool has_checkpoint_data     = m_avg_type!=OutputAvgType::Instant && not output_every_step;
  const bool is_full_checkpoint_step = is_checkpoint_step && has_checkpoint_data && not is_output_step;
  const bool is_write_step           = is_output_step || is_checkpoint_step;

  // Update counters
  ++m_output_control.nsamples_since_last_write;
  if (not is_t0_output) {
    // In case REST_OPT=nsteps, don't count t0 output as one of those steps
    // NOTE: for m_output_control, it doesn't matter, since it'll be reset to 0 before we return
    ++m_checkpoint_control.nsamples_since_last_write;
  }

  // Create and setup output/checkpoint file(s), if necessary
  start_timer(timer_root+"::get_new_file");
  auto setup_output_file = [&](IOControl& control, IOFileSpecs& filespecs) {
    // Check if the new snapshot fits, if not, close the file
    // NOTE: if output is average/max/min AND we save one file per month/year,
    //       we don't want to check if *this* timestamp fits in the file, since
    //       it may be in the next month/year even though most of the time averaging
    //       window is in the right year. E.g., if the avg is over the whole month,
    //       you would end up saving Jan average in the Feb file, since the end
    //       of the window is Feb 1st 00:00:00. So instead, we use the *start*
    //       of the avg window. If the avg is such that it spans 2 months (e.g.,
    //       a 7-day avg), then where we put the avg is arbitrary, so our choice
    //       is still fine.
    util::TimeStamp snapshot_start;
    if (m_avg_type==OutputAvgType::Instant or filespecs.storage.type==NumSnaps) {
      snapshot_start = timestamp;
    } else {
      snapshot_start = m_case_t0;
      snapshot_start += m_time_bnds[0];
    }

    // Check if we need to open a new file
    if (not filespecs.is_open) {
      filespecs.filename = compute_filename (filespecs,timestamp);
      // Register all dims/vars, write geometry data (e.g. lat/lon/hyam/hybm)
      setup_file(filespecs,control);
    }

    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (m_io_comm.am_i_root() and filespecs.is_restart_file()) {
      std::ofstream rpointer;
      if (m_is_model_restart_output) {
        rpointer.open("rpointer.atm");  // Open rpointer and nuke its content
      } else if (is_checkpoint_step) {
        // Output restart unit tests do not have a model-output stream that generates rpointer.atm,
        // so allow to skip the next check for them.
        auto is_unit_testing = m_params.sublist("Checkpoint Control").get("is_unit_testing",false);
        EKAT_REQUIRE_MSG (is_unit_testing || std::ifstream("rpointer.atm").good(),
            "Error! Cannot find rpointer.atm file to append history restart file in.\n"
            " Model restart output is supposed to be in charge of creating rpointer.atm.\n"
            " There are two possible causes:\n"
            "   1. You have a 'Checkpoint Control' list in your output stream, but no Scorpio::model_restart\n"
            "      section in the input yaml file. This makes no sense, please correct.\n"
            "   2. The current implementation assumes that the model restart OutputManager runs\n"
            "      *before* any other output stream (so it can nuke rpointer.atm if already existing).\n"
            "      If this has changed, we need to revisit this piece of the code.\n");
        rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append to it
      }
      rpointer << filespecs.filename << std::endl;
    }

    if (m_atm_logger) {
      m_atm_logger->info("[EAMxx::output_manager] - Writing " + e2str(filespecs.ftype) + ":");
      m_atm_logger->info("[EAMxx::output_manager]      FILE: " + filespecs.filename);
    }
  };

  if (is_output_step) {
    setup_output_file(m_output_control,m_output_file_specs);

    // Update time (must be done _before_ writing fields)
    update_time(m_output_file_specs.filename,timestamp.days_from(m_case_t0));
  }
  if (is_checkpoint_step) {
    setup_output_file(m_checkpoint_control,m_checkpoint_file_specs);

    if (is_full_checkpoint_step) {
      // Update time (must be done _before_ writing fields)
      update_time(m_checkpoint_file_specs.filename,timestamp.days_from(m_case_t0));
    }
  }
  stop_timer(timer_root+"::get_new_file");

  // Run the output streams
  start_timer(timer_root+"::run_output_streams");
  const auto& fields_write_filename = is_output_step ? m_output_file_specs.filename : m_checkpoint_file_specs.filename;
  for (auto& it : m_output_streams) {
    // Note: filename only matters if is_output_step || is_full_checkpoint_step=true. In that case, it will definitely point to a valid file name.
    if (m_atm_logger) {
      m_atm_logger->debug("[OutputManager]: writing fields from grid " + it->get_io_grid()->name() + "...\n");
    }
    it->run(fields_write_filename,is_output_step,is_full_checkpoint_step,m_output_control.nsamples_since_last_write,is_t0_output);
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
      if (m_atm_logger) {
        m_atm_logger->debug("[OutputManager]: writing globals...\n");
      }

      // Since we wrote to file we need to reset the timestamps
      control.last_write_ts = timestamp;
      control.compute_next_write_ts();
      control.nsamples_since_last_write = 0;

      if (m_is_model_restart_output) {
        // Only write nsteps on model restart
        set_attribute(filespecs.filename,"GLOBAL","nsteps",timestamp.get_num_steps());
      } else {
        if (filespecs.ftype==FileType::HistoryRestart) {
          // Update the date of last write and sample size
          write_timestamp (filespecs.filename,"last_write",m_output_control.last_write_ts,true);
          scorpio::set_attribute (filespecs.filename,"GLOBAL","last_output_filename",m_output_file_specs.filename);
          scorpio::set_attribute (filespecs.filename,"GLOBAL","num_snapshots_since_last_write",m_output_control.nsamples_since_last_write);
          scorpio::set_attribute (filespecs.filename,"GLOBAL","last_output_file_num_snaps",m_output_file_specs.storage.num_snapshots_in_file);
        }
        // Write these in both output and rhist file. The former, b/c we need these info when we postprocess
        // output, and the latter b/c we want to make sure these params don't change across restarts
        set_attribute(filespecs.filename,"GLOBAL","averaging_type",e2str(m_avg_type));
        set_attribute(filespecs.filename,"GLOBAL","averaging_frequency_units",m_output_control.frequency_units);
        set_attribute(filespecs.filename,"GLOBAL","averaging_frequency",m_output_control.frequency);
        set_attribute(filespecs.filename,"GLOBAL","file_max_storage_type",e2str(m_output_file_specs.storage.type));
        if (m_output_file_specs.storage.type==NumSnaps) {
          set_attribute(filespecs.filename,"GLOBAL","max_snapshots_per_file",m_output_file_specs.storage.max_snapshots_in_file);
        }
        const auto& fp_precision = m_params.get<std::string>("Floating Point Precision");
        set_attribute(filespecs.filename,"GLOBAL","fp_precision",fp_precision);
      }

      // Write all stored globals
      for (const auto& it : m_globals) {
        const auto& name = it.first;
        const auto& any = it.second;
        if (any.isType<int>()) {
          set_attribute(filespecs.filename,"GLOBAL",name,ekat::any_cast<int>(any));
        } else if (any.isType<std::int64_t>()) {
          set_attribute(filespecs.filename,"GLOBAL",name,ekat::any_cast<std::int64_t>(any));
        } else if (any.isType<float>()) {
          set_attribute(filespecs.filename,"GLOBAL",name,ekat::any_cast<float>(any));
        } else if (any.isType<double>()) {
          set_attribute(filespecs.filename,"GLOBAL",name,ekat::any_cast<double>(any));
        } else if (any.isType<std::string>()) {
          set_attribute(filespecs.filename,"GLOBAL",name,ekat::any_cast<std::string>(any));
        } else {
          EKAT_ERROR_MSG (
              "Error! Invalid concrete type for IO global.\n"
              " - global name: " + it.first + "\n"
              " - type id    : " + any.content().type().name() + "\n");
        }
      }

      // We're adding one snapshot to the file
      filespecs.storage.update_storage(timestamp);

      // NOTE: for checkpoint files, unless we write restart data, we did not update time,
      //       which means we cannot write any variable (the check var.num_records==time.length
      //       would fail)
      if (m_time_bnds.size()>0 and
          (filespecs.ftype!=FileType::HistoryRestart or is_full_checkpoint_step)) {
        scorpio::write_var(filespecs.filename, "time_bnds", m_time_bnds.data());
      }

      close_or_flush_if_needed(filespecs,control);
    };

    start_timer(timer_root+"::update_snapshot_tally");
    // Important! Process output file first, and hist restart (if any) second.
    // That's b/c write_global_data will update m_output_control.last_write_ts,
    // which is later written as global data in the hist restart file
    if (is_output_step) {
      write_global_data(m_output_control,m_output_file_specs);
    }
    if (is_checkpoint_step) {
      write_global_data(m_checkpoint_control,m_checkpoint_file_specs);

      // Always flush output during checkpoints (assuming we opened it already)
      if (m_output_file_specs.is_open) {
        scorpio::flush_file (m_output_file_specs.filename);
      }
    }
    stop_timer(timer_root+"::update_snapshot_tally");
    if (is_output_step && m_time_bnds.size()>0) {
      m_time_bnds[0] = m_time_bnds[1];
    }
  }

  stop_timer("EAMxx::IO::" + m_params.name());
  stop_timer(timer_root);
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  // Close any output file still open
  if (m_output_file_specs.is_open) {
    scorpio::release_file (m_output_file_specs.filename);
  }
  if (m_checkpoint_file_specs.is_open) {
    scorpio::release_file (m_checkpoint_file_specs.filename);
  }

  // Reset everything to a default constructed object.
  // NOTE: it's themptying to std::swap(*this,OutputManager()),
  //       but that calls ~OutputManager() on the destructor,
  //       which in turns calls finalize, causing endless recursion.
  m_output_streams = {};
  m_geo_data_streams = {};
  m_globals.clear();
  m_io_comm = {};
  m_params  = {};
  m_filename_prefix = {};
  m_time_bnds = {};
  m_avg_type = {};
  m_output_control = {};
  m_checkpoint_control = {};
  m_output_file_specs = {};
  m_checkpoint_file_specs = {};
  m_case_t0 = {};
  m_run_t0 = {};
  m_atm_logger = {};
}

long long OutputManager::res_dep_memory_footprint () const {
  long long mf = 0;
  for (const auto& os : m_output_streams) {
    mf += os->res_dep_memory_footprint();
  }

  return mf;
}

std::string OutputManager::
compute_filename (const IOFileSpecs& file_specs,
                  const util::TimeStamp& timestamp) const
{
  auto filename = m_filename_prefix + file_specs.suffix();
  const auto& control = m_output_control;

  // Always add avg type and frequency info
  filename += "." + e2str(m_avg_type);
  filename += "." + control.frequency_units+ "_x" + std::to_string(control.frequency);

  // For standalone EAMxx, we may have 2+ versions of the same test running with two
  // different choices of ranks. To avoid name clashing for the output files,
  // add the comm size to the output file name.
  if (is_scream_standalone()) {
    filename += ".np" + std::to_string(m_io_comm.size());
  }

  // Always add a time stamp
  auto ts = (m_avg_type==OutputAvgType::Instant || file_specs.ftype==FileType::HistoryRestart)
          ? timestamp : control.last_write_ts;

  int ts_string_len = 0;
  switch (file_specs.storage.type) {
    case Yearly:   ts_string_len = 4;  break; // YYYY
    case Monthly:  ts_string_len = 7;  break; // YYYY-MM
    case Daily:    ts_string_len = 10; break; // YYYY-MM-DD
    case NumSnaps: ts_string_len = 16; break; // YYYY-MM-DD-XXXXX
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
  }
  filename += "." + ts.to_string().substr(0,ts_string_len);

  return filename + ".nc";
}

void OutputManager::
setup_internals (const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs)
{
  using vos_t = std::vector<std::string>;

  if (m_is_model_restart_output) {
    // We build some restart parameters internally
    m_avg_type = OutputAvgType::Instant;
    m_output_file_specs.storage.type = NumSnaps;
    m_output_file_specs.storage.max_snapshots_in_file = 1;
    m_output_file_specs.flush_frequency = 1;

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
    m_filename_prefix = m_params.get<std::string>("filename_prefix");

    // Hard code some parameters in case we access them later
    m_params.set<std::string>("Floating Point Precision","real");
  } else {
    auto avg_type = m_params.get<std::string>("Averaging Type");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
        "Error! Unsupported averaging type '" + avg_type + "'.\n"
        "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

    const auto& storage_type = m_params.get<std::string>("file_max_storage_type","num_snapshots");
    auto& storage = m_output_file_specs.storage;
    constexpr auto large_int = std::numeric_limits<int>::max();
    if (storage_type=="num_snapshots") {
      storage.type = NumSnaps;
      storage.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File",large_int);
      EKAT_REQUIRE_MSG (storage.max_snapshots_in_file>0,
          "Error! Value for 'Max Snapshots Per File' should be positive.\n"
          "       To request 'unlimited' storage, leave parameter unset.\n");
    } else if (storage_type=="one_year") {
      storage.type = Yearly;
    } else if (storage_type=="one_month") {
      storage.type = Monthly;
    } else if (storage_type=="one_day") {
      storage.type = Daily;
    } else {
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
    }
    m_filename_prefix = m_params.get<std::string>("filename_prefix");
    m_output_file_specs.flush_frequency = m_params.get("flush_frequency",large_int);

    // Allow user to ask for higher precision for normal model output,
    // but default to single to save on storage
    const auto& prec = m_params.get<std::string>("Floating Point Precision", "single");
    vos_t valid_prec = {"single", "float", "double", "real"};
    EKAT_REQUIRE_MSG (ekat::contains(valid_prec,prec),
        "Error! Invalid/unsupported value for 'Floating Point Precision'.\n"
        "  - input value: " + prec + "\n"
        "  - supported values: float, single, double, real\n");
  }

  // Output control
  EKAT_REQUIRE_MSG(m_params.isSublist("output_control"),
      "Error! The output control YAML file for " + m_filename_prefix + " is missing the sublist 'output_control'");
  auto& out_control_pl = m_params.sublist("output_control");
  m_output_control.set_frequency_units(out_control_pl.get<std::string>("frequency_units"));

  // In case output is disabled, no point in doing anything else
  if (m_output_control.frequency_units=="none" || m_output_control.frequency_units=="never") {
    return;
  }
  m_output_control.frequency = out_control_pl.get<int>("Frequency");
  EKAT_REQUIRE_MSG (m_output_control.frequency>0,
      "Error! Invalid frequency (" + std::to_string(m_output_control.frequency) + ") in Output Control. Please, use positive number.\n");

  // Determine which timestamp to use a reference for output frequency.  Two options:
  // 	1. use_case_as_start_reference: TRUE  - implies we want to calculate frequency from the beginning of the whole simulation, even if this is a restarted run.
  // 	2. use_case_as_start_reference: FALSE - implies we want to base the frequency of output on when this particular simulation started.
  // Note, (2) is needed for restarts since the restart frequency in CIME assumes a reference of when this run began.
  const bool use_case_as_ref = out_control_pl.get<bool>("use_case_as_start_reference",!m_is_model_restart_output);
  const bool perform_history_restart = m_params.sublist("Restart").get("Perform Restart",true);
  const auto& start_ref = use_case_as_ref and perform_history_restart ? m_case_t0 : m_run_t0;
  m_output_control.last_write_ts = start_ref;
  m_output_control.compute_next_write_ts();

  // File specs
  m_save_grid_data = out_control_pl.get("save_grid_data",!m_is_model_restart_output);
  m_output_file_specs.ftype = m_is_model_restart_output ? FileType::ModelRestart : FileType::ModelOutput;

  if (m_params.isSublist("Checkpoint Control")) {
    auto& pl = m_params.sublist("Checkpoint Control");
    m_checkpoint_control.set_frequency_units(pl.get<std::string>("frequency_units"));

    if (m_checkpoint_control.output_enabled()) {
      m_checkpoint_control.frequency = pl.get<int>("Frequency");
      EKAT_REQUIRE_MSG (m_output_control.frequency>0,
          "Error! Invalid frequency (" + std::to_string(m_checkpoint_control.frequency) + ") in Checkpoint Control. Please, use positive number.\n");

      m_checkpoint_control.last_write_ts = m_run_t0;
      m_checkpoint_control.compute_next_write_ts();

      // File specs
      m_checkpoint_file_specs.storage.type = NumSnaps;
      m_checkpoint_file_specs.storage.max_snapshots_in_file = 1;
      m_checkpoint_file_specs.flush_frequency = 1;
      m_checkpoint_file_specs.ftype = FileType::HistoryRestart;
    }
  }

  // Set the iotype to use for the output file
  std::string iotype = m_params.get<std::string>("iotype", "default");
  m_output_file_specs.iotype = scorpio::str2iotype(iotype);
  m_checkpoint_file_specs.iotype = scorpio::str2iotype(iotype);
}

/*===============================================================================================*/
void OutputManager::
setup_file (      IOFileSpecs& filespecs,
            const IOControl& control)
{
  const bool is_checkpoint_step = &control==&m_checkpoint_control;

  std::string fp_precision = is_checkpoint_step
                           ? "real"
                           : m_params.get<std::string>("Floating Point Precision");

  const auto& filename = filespecs.filename;
  // Register new netCDF file for output. Check if we need to append to an existing file
  auto mode = m_resume_output_file ? scorpio::Append : scorpio::Write;
  scorpio::register_file(filename,mode,filespecs.iotype);
  if (m_resume_output_file) {
    // We may have resumed an output file that contains extra snapshots *after* the restart time.
    // E.g., if we output every step and the run crashed a few steps after writing the restart.
    // In that case, we need to reset the time dimension in the output file, so that the extra
    // snapshots will be overwritten.
    const auto all_times = scorpio::get_all_times(filename);
    int ntimes = all_times.size();
    int ngood  = 0;
    for (const auto& t : all_times) {
      auto keep = t<=m_output_control.last_write_ts.days_from(m_case_t0);
      if (keep) {
        ++ngood;
      } else {
        break;
      }
    }
    if (ngood<ntimes) {
      scorpio::reset_unlimited_dim_len(filename,ngood);
    }
    scorpio::redef(filename);
  } else {
    // Register time (and possibly time_bnds) var(s)
    auto time_units="days since " + m_case_t0.get_date_string() + " " + m_case_t0.get_time_string();
    scorpio::define_time(filename,time_units,"time");

    scorpio::define_var(filename,"time",time_units,{}, "double", "double",true);
    if (use_leap_year()) {
      scorpio::set_attribute (filename,"time","calendar","gregorian");
    } else {
      scorpio::set_attribute (filename,"time","calendar","noleap");
    }

    if (m_avg_type!=OutputAvgType::Instant) {
      // First, ensure a 'dim2' dimension with len=2 is registered.
      scorpio::define_dim(filename,"dim2",2);
      scorpio::define_var(filename,"time_bnds",time_units,{"dim2"},"double","double",true);

      // Make it clear how the time_bnds should be interpreted
      scorpio::set_attribute<std::string> (filename,"time_bnds","note","right endpoint accumulation");

      // I'm not sure what's the point of this, but CF conventions seem to require it
      scorpio::set_attribute<std::string> (filename,"time","bounds","time_bnds");
    }

    write_timestamp(filename,"case_t0",m_case_t0);
    write_timestamp(filename,"run_t0",m_run_t0);
    scorpio::set_attribute(filename,"GLOBAL","averaging_type",e2str(m_avg_type));
    scorpio::set_attribute(filename,"GLOBAL","averaging_frequency_units",m_output_control.frequency_units);
    scorpio::set_attribute(filename,"GLOBAL","averaging_frequency",m_output_control.frequency);
    scorpio::set_attribute(filespecs.filename,"GLOBAL","file_max_storage_type",e2str(m_output_file_specs.storage.type));
    if (m_output_file_specs.storage.type==NumSnaps) {
      scorpio::set_attribute(filename,"GLOBAL","max_snapshots_per_file",m_output_file_specs.storage.max_snapshots_in_file);
    }
    scorpio::set_attribute(filename,"GLOBAL","fp_precision",fp_precision);
    set_file_header(filespecs);
  }

  // Make all output streams register their dims/vars
  for (auto& it : m_output_streams) {
    it->setup_output_file(filename,fp_precision,mode);
  }

  // If grid data is needed,  also register geo data fields. Skip if file is resumed,
  // since grid data was written in the previous run
  if (m_save_grid_data and not filespecs.is_restart_file() and not m_resume_output_file) {
    for (auto& it : m_geo_data_streams) {
      it->setup_output_file(filename,fp_precision,mode);
    }
  }

  scorpio::enddef (filename);

  if (m_save_grid_data and not filespecs.is_restart_file() and not m_resume_output_file) {
    // Immediately run the geo data streams
    for (const auto& it : m_geo_data_streams) {
      it->run(filename,true,false,0);
    }
  }

  filespecs.is_open = true;

  m_resume_output_file = false;
}
/*===============================================================================================*/
void OutputManager::set_file_header(const IOFileSpecs& file_specs)
{
  auto& p = m_params.sublist("provenance");
  auto now = std::chrono::system_clock::now();
  std::time_t time = std::chrono::system_clock::to_time_t(now);
  std::stringstream timestamp;
  timestamp << "created on " << std::ctime(&time);
  std::string ts_str = timestamp.str();
  ts_str = std::strtok(&ts_str[0],"\n"); // Remove the \n appended by ctime

  const auto& filename = file_specs.filename;

  auto set_str_att = [&](const std::string& name, const std::string& val) {
    scorpio::set_attribute(filename,"GLOBAL",name,val);
  };

  set_str_att("case",p.get<std::string>("caseid","NONE"));
  set_str_att("source","E3SM Atmosphere Model (EAMxx)");
  set_str_att("eamxx_version",EAMXX_VERSION);
  set_str_att("git_version",p.get<std::string>("git_version",EAMXX_GIT_VERSION));
  set_str_att("hostname",p.get<std::string>("hostname","UNKNOWN"));
  set_str_att("username",p.get<std::string>("username","UNKNOWN"));
  set_str_att("atm_initial_conditions_file",p.get<std::string>("initial_conditions_file","NONE"));
  set_str_att("topography_file",p.get<std::string>("topography_file","NONE"));
  set_str_att("contact","e3sm-data-support@llnl.gov");
  set_str_att("institution_id","E3SM-Project");
  set_str_att("realm","atmos");
  set_str_att("history",ts_str);
  set_str_att("Conventions","CF-1.8");
  set_str_att("product",e2str(file_specs.ftype));
}
void OutputManager::
close_or_flush_if_needed (      IOFileSpecs& file_specs,
                          const IOControl&   control) const
{
  if (not file_specs.storage.snapshot_fits(control.next_write_ts)) {
    scorpio::release_file(file_specs.filename);
    file_specs.close();
  } else if (file_specs.file_needs_flush()) {
    scorpio::flush_file (file_specs.filename);
  }
}

void OutputManager::
push_to_logger()
{
  // If no atm logger set then don't do anything
  if (!m_atm_logger) return;

  auto bool_to_string = [](const bool x) {
    std::string y = x ? "YES" : "NO";
    return y;
  };

  auto rt_to_string = [](RunType rt) {
    std::string s = rt==RunType::Initial ? "Initial" : "Restart";
    return s;
  };

  m_atm_logger->info("[EAMxx::output_manager] - New Output stream");
  m_atm_logger->info("           Filename prefix: " + m_filename_prefix);
  m_atm_logger->info("                    Run t0: " + m_run_t0.to_string());
  m_atm_logger->info("                   Case t0: " + m_case_t0.to_string());
  m_atm_logger->info("              Reference t0: " + m_output_control.last_write_ts.to_string());
  m_atm_logger->info("         Is Restart File ?: " + bool_to_string(m_is_model_restart_output));
  m_atm_logger->info("                 Run type : " + rt_to_string(m_run_type));
  m_atm_logger->info("            Averaging Type: " + e2str(m_avg_type));
  m_atm_logger->info("          Output Frequency: " + std::to_string(m_output_control.frequency) + " " + m_output_control.frequency_units);
  switch (m_output_file_specs.storage.type) {
    case NumSnaps:
    {
      auto ms = m_output_file_specs.storage.max_snapshots_in_file;
      m_atm_logger->info("             File Capacity: " + (ms>0 ? std::to_string(ms) + "snapshots" : "UNLIMITED"));
      break;
    }
    case Monthly:
      m_atm_logger->info("             File Capacity: one month per file");
      break;
    case Yearly:
      m_atm_logger->info("             File Capacity: one year per file");
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
  }
  m_atm_logger->info("      Includes Grid Data ?: " + bool_to_string(m_save_grid_data));
  // List each GRID - TODO
  // List all FIELDS - TODO
}

} // namespace scream
