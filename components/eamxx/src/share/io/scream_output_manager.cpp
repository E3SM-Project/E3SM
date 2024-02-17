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

  // Read input parameters and setup internal data
  set_params(params,field_mgrs);

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
  if (m_is_restarted_run and not m_is_model_restart_output) {
    // Allow to skip history restart, or to specify a filename_prefix for the restart file
    // that is different from the filename_prefix of the current output.
    auto& restart_pl = m_params.sublist("Restart");
    bool perform_history_restart = restart_pl.get("Perform Restart",true);
    auto hist_restart_filename_prefix = restart_pl.get("filename_prefix",m_filename_prefix);

    if (perform_history_restart) {
      using namespace scorpio;
      auto rhist_file = find_filename_in_rpointer(hist_restart_filename_prefix,false,m_io_comm,m_run_t0);

      // From restart file, get the time of last write, as well as the current size of the avg sample
      m_output_control.last_write_ts = read_timestamp(rhist_file,"last_write");
      m_output_control.compute_next_write_ts();
      m_output_control.nsamples_since_last_write = get_attribute<int>(rhist_file,"num_snapshots_since_last_write");

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
      auto old_freq = scorpio::get_attribute<int>(rhist_file,"averaging_frequency");
      EKAT_REQUIRE_MSG (old_freq == m_output_control.frequency,
          "Error! Cannot change frequency when performing history restart.\n"
          "  - old freq: " << old_freq << "\n"
          "  - new freq: " << m_output_control.frequency << "\n");
      auto old_freq_units = scorpio::get_attribute<std::string>(rhist_file,"averaging_frequency_units");
      EKAT_REQUIRE_MSG (old_freq_units == m_output_control.frequency_units,
          "Error! Cannot change frequency units when performing history restart.\n"
          "  - old freq units: " << old_freq_units << "\n"
          "  - new freq units: " << m_output_control.frequency_units << "\n");
      auto old_avg_type = scorpio::get_attribute<std::string>(rhist_file,"averaging_type");
      EKAT_REQUIRE_MSG (old_avg_type == e2str(m_avg_type),
          "Error! Cannot change avg type when performing history restart.\n"
          "  - old avg type: " << old_avg_type + "\n"
          "  - new avg type: " << e2str(m_avg_type) << "\n");


      auto old_storage_type = scorpio::get_attribute<std::string>(rhist_file,"file_max_storage_type");
      EKAT_REQUIRE_MSG (old_storage_type == e2str(m_output_file_specs.storage.type),
          "Error! Cannot change file storage type when performing history restart.\n"
          "  - old file_max_storage_type: " << old_storage_type << "\n"
          "  - new file_max_storage_type: " << e2str(m_output_file_specs.storage.type) << "\n"
          "If you *really* want to change the file storage type, you need to force using a new file, setting\n"
          "  Restart:\n"
          "    force_new_file: true\n");
      if (old_storage_type=="num_snapshot") {
        auto old_max_snaps = scorpio::get_attribute<int>(rhist_file,"max_snapshots_per_file");
        EKAT_REQUIRE_MSG (old_max_snaps == m_output_file_specs.storage.max_snapshots_in_file,
            "Error! Cannot change max snapshots per file when performing history restart.\n"
            "  - old max snaps: " << old_max_snaps << "\n"
            "  - new max snaps: " << m_output_file_specs.storage.max_snapshots_in_file << "\n"
            "If you *really* want to change the file capacity, you need to force using a new file, setting\n"
            "  Restart:\n"
            "    force_new_file: true\n");
      }
      std::string fp_precision = m_params.get<std::string>("Floating Point Precision");
      auto old_fp_precision = scorpio::get_attribute<std::string>(rhist_file,"fp_precision");
      EKAT_REQUIRE_MSG (old_fp_precision == fp_precision,
          "Error! Cannot change floating point precision when performing history restart.\n"
          "  - old fp precision: " << old_fp_precision << "\n"
          "  - new fp precision: " << fp_precision << "\n");

      // Check if the prev run wrote any output file (it may have not, if the restart was written
      // before the 1st output step). If there is a file, check if there's still room in it.
      const auto& last_output_filename = get_attribute<std::string>(rhist_file,"last_output_filename");
      m_resume_output_file = last_output_filename!="" and not restart_pl.get("force_new_file",false);
      if (m_resume_output_file) {

        scorpio::register_file(last_output_filename,scorpio::Read);
        int num_snaps = scorpio::get_dimlen(last_output_filename,"time");
        scorpio::eam_pio_closefile(last_output_filename);

        m_output_file_specs.filename = last_output_filename;
        m_output_file_specs.is_open = true;
        m_output_file_specs.storage.num_snapshots_in_file = num_snaps;
        // The setup_file call will not register any new variable (the file is in Append mode,
        // so all dims/vars must already be in the file). However, it will register decompositions,
        // since those are a property of the run, not of the file.
        setup_file(m_output_file_specs,m_output_control);
      }
    }
  }

  // If m_time_bnds.size()>0, it was already inited during restart
  if (m_avg_type!=OutputAvgType::Instant && m_time_bnds.size()==0) {
    // Init the left hand point of time_bnds based on run/case t0.
    m_time_bnds.resize(2);
    m_time_bnds[0] = m_run_t0.days_from(m_case_t0);
  } else if (m_output_control.output_enabled() && m_run_t0==m_case_t0 && !m_is_model_restart_output) {
    // In order to trigger a t0 write, we need to have next_write_ts matching run_t0
    m_output_control.next_write_ts = m_run_t0;
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
  if (not m_output_control.output_enabled()) {
    return;
  }

  if (m_atm_logger) {
    m_atm_logger->debug("[OutputManager::run] filename_prefix: " + m_filename_prefix + "\n");
  }

  using namespace scorpio;

  // We did not have dt at init time. Now we can compute it based on timestamp and the ts of last write.
  // NOTE: dt is only needed for frequency_units='nsteps'. For other units, the function returns immediately
  m_output_control.compute_dt(timestamp);
  if (m_checkpoint_control.output_enabled()) {
    m_checkpoint_control.compute_dt(timestamp);
  }

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
    if (not filespecs.storage.snapshot_fits(snapshot_start)) {
      eam_pio_closefile(filespecs.filename);
      filespecs.close();
    }

    // Check if we need to open a new file
    if (not filespecs.is_open) {
      filespecs.filename = compute_filename (control,filespecs,timestamp);
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
    pio_update_time(m_output_file_specs.filename,timestamp.days_from(m_case_t0));
  }
  if (is_checkpoint_step) {
    setup_output_file(m_checkpoint_control,m_checkpoint_file_specs);

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
      control.update_write_timestamps();

      if (m_is_model_restart_output) {
        // Only write nsteps on model restart
        set_attribute(filespecs.filename,"nsteps",timestamp.get_num_steps());
      } else {
        if (filespecs.ftype==FileType::HistoryRestart) {
          // Update the date of last write and sample size
          scorpio::write_timestamp (filespecs.filename,"last_write",m_output_control.last_write_ts);
          scorpio::set_attribute (filespecs.filename,"last_write_nsteps",m_output_control.last_write_ts.get_num_steps());
          scorpio::set_attribute (filespecs.filename,"last_output_filename",m_output_file_specs.filename);
          scorpio::set_attribute (filespecs.filename,"num_snapshots_since_last_write",m_output_control.nsamples_since_last_write);
        }
        // Write these in both output and rhist file. The former, b/c we need these info when we postprocess
        // output, and the latter b/c we want to make sure these params don't change across restarts
        set_attribute(filespecs.filename,"averaging_type",e2str(m_avg_type));
        set_attribute(filespecs.filename,"averaging_frequency_units",m_output_control.frequency_units);
        set_attribute(filespecs.filename,"averaging_frequency",m_output_control.frequency);
        set_attribute(filespecs.filename,"file_max_storage_type",e2str(m_output_file_specs.storage.type));
        if (m_output_file_specs.storage.type==NumSnaps) {
          set_attribute(filespecs.filename,"max_snapshots_per_file",m_output_file_specs.storage.max_snapshots_in_file);
        }
        const auto& fp_precision = m_params.get<std::string>("Floating Point Precision");
        set_attribute(filespecs.filename,"fp_precision",fp_precision);
      }

      // Write all stored globals
      for (const auto& it : m_globals) {
        const auto& name = it.first;
        const auto& any = it.second;
        set_any_attribute(filespecs.filename,name,any);
      }

      // We're adding one snapshot to the file
      filespecs.storage.update_storage(timestamp);

      if (m_time_bnds.size()>0) {
        scorpio::grid_write_data_array(filespecs.filename, "time_bnds", m_time_bnds.data(), 2);
      }

      // Check if we need to flush the output file
      if (filespecs.file_needs_flush()) {
        eam_flush_file (filespecs.filename);
      }
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
    }
    stop_timer(timer_root+"::update_snapshot_tally");
    if (is_output_step && m_time_bnds.size()>0) {
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
                  const util::TimeStamp& timestamp) const
{
  auto filename = m_filename_prefix + file_specs.suffix();

  // Always add avg type and frequency info
  filename += "." + e2str(m_avg_type);
  filename += "." + control.frequency_units+ "_x" + std::to_string(control.frequency);

  // Optionally, add number of mpi ranks (useful mostly in unit tests, to run multiple MPI configs in parallel)
  if (m_params.get<bool>("MPI Ranks in Filename")) {
    filename += ".np" + std::to_string(m_io_comm.size());
  }

  // Always add a time stamp
  auto ts = (m_avg_type==OutputAvgType::Instant || file_specs.ftype==FileType::HistoryRestart)
          ? timestamp : control.last_write_ts;

  switch (file_specs.storage.type) {
    case NumSnaps:
      filename += "." + ts.to_string(); break;
    case Yearly:
      filename += "." + std::to_string(ts.get_year()); break;
    case Monthly:
      filename += "." + std::to_string(ts.get_year()) + "-" + std::to_string(ts.get_month()); break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported file storage type.\n");
  }

  return filename + ".nc";
}

void OutputManager::
set_params (const ekat::ParameterList& params,
            const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs)
{
  using vos_t = std::vector<std::string>;

  m_params = params;

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
    m_params.set("MPI Ranks in Filename",false);
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
    } else if (storage_type=="one_year") {
      storage.type = Yearly;
    } else if (storage_type=="one_month") {
      storage.type = Monthly;
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

    // If not set, hard code to false for CIME cases, and true for standalone,
    // since standalone may be running multiple versions of the same test at once
    if (not m_params.isParameter("MPI Ranks in Filename")) {
      m_params.set("MPI Ranks in Filename",is_scream_standalone());
    }
  }

  // Output control
  EKAT_REQUIRE_MSG(m_params.isSublist("output_control"),
      "Error! The output control YAML file for " + m_filename_prefix + " is missing the sublist 'output_control'");
  auto& out_control_pl = m_params.sublist("output_control");
  m_output_control.frequency_units = out_control_pl.get<std::string>("frequency_units");

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
  // NOTE: m_is_restarted_run and not m_is_model_restart_output, these timestamps will be corrected once we open the hist restart file
  const bool use_case_as_ref = out_control_pl.get<bool>("use_case_as_start_reference",!m_is_model_restart_output);
  const auto& start_ref = use_case_as_ref ? m_case_t0 : m_run_t0;
  m_output_control.last_write_ts = start_ref;
  m_output_control.compute_next_write_ts();

  // File specs
  m_save_grid_data = out_control_pl.get("save_grid_data",!m_is_model_restart_output);
  m_output_file_specs.ftype = m_is_model_restart_output ? FileType::ModelRestart : FileType::ModelOutput;

  if (m_params.isSublist("Checkpoint Control")) {
    auto& pl = m_params.sublist("Checkpoint Control");
    m_checkpoint_control.frequency_units = pl.get<std::string>("frequency_units");

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

}

/*===============================================================================================*/
void OutputManager::
setup_file (      IOFileSpecs& filespecs,
            const IOControl& control)
{
  using namespace scorpio;

  const bool is_checkpoint_step = &control==&m_checkpoint_control;

  std::string fp_precision = is_checkpoint_step
                           ? "real"
                           : m_params.get<std::string>("Floating Point Precision");

  const auto& filename = filespecs.filename;
  // Register new netCDF file for output. Check if we need to append to an existing file
  auto mode = m_resume_output_file ? Append : Write;
  register_file(filename,mode);
  if (m_resume_output_file) {
    eam_pio_redef(filename);
  }

  // Note: length=0 is how scorpio recognizes that this is an 'unlimited' dimension, which
  // allows to write as many timesnaps as we desire.
  register_dimension(filename,"time","time",0,false);

  // Register time (and possibly time_bnds) var(s)
  auto time_units="days since " + m_case_t0.get_date_string() + " " + m_case_t0.get_time_string();
  register_variable(filename,"time","time",time_units,{"time"}, "double", "double","time");
  if (use_leap_year()) {
    set_variable_metadata (filename,"time","calendar","gregorian");
  } else {
    set_variable_metadata (filename,"time","calendar","noleap");
  }
  if (m_avg_type!=OutputAvgType::Instant) {
    // First, ensure a 'dim2' dimension with len=2 is registered.
    register_dimension(filename,"dim2","dim2",2,false);
    register_variable(filename,"time_bnds","time_bnds",time_units,{"dim2","time"},"double","double","time-dim2");
    
    // Make it clear how the time_bnds should be interpreted
    set_variable_metadata(filename,"time_bnds","note","right endpoint accumulation");

    // I'm not sure what's the point of this, but CF conventions seem to require it
    set_variable_metadata (filename,"time","bounds","time_bnds");
  }

  if (not m_resume_output_file) {
    // Finish the definition phase for this file.
    write_timestamp(filename,"case_t0",m_case_t0);
    write_timestamp(filename,"run_t0",m_run_t0);
    set_attribute(filename,"averaging_type",e2str(m_avg_type));
    set_attribute(filename,"averaging_frequency_units",m_output_control.frequency_units);
    set_attribute(filename,"averaging_frequency",m_output_control.frequency);
    set_attribute(filespecs.filename,"file_max_storage_type",e2str(m_output_file_specs.storage.type));
    if (m_output_file_specs.storage.type==NumSnaps) {
      set_attribute(filename,"max_snapshots_per_file",m_output_file_specs.storage.max_snapshots_in_file);
    }
    set_attribute(filename,"fp_precision",fp_precision);
    set_file_header(filespecs);
  }

  // Set degree of freedom for "time" and "time_bnds"
  scorpio::offset_t time_dof[1] = {0};
  set_dof(filename,"time",1,time_dof);
  if (m_avg_type!=OutputAvgType::Instant) {
    scorpio::offset_t time_bnds_dofs[2] = {0,1};
    set_dof(filename,"time_bnds",2,time_bnds_dofs);
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

  // When resuming a file, PIO opens it in data mode.
  // NOTE: all the above register_dimension/register_variable are already checking that
  //       the dims/vars are already in the file (we don't allow adding dims/vars)
  eam_pio_enddef (filename);

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
  using namespace scorpio;

  // TODO: All attributes marked TODO below need to be set.  Hopefully by a universal value that reflects
  // what the attribute is.  For example, git-hash should be the git-hash associated with this version of
  // the code at build time for this executable.
  auto& p = m_params.sublist("provenance");
  auto now = std::chrono::system_clock::now();
  std::time_t time = std::chrono::system_clock::to_time_t(now);
  std::stringstream timestamp;
  timestamp << "created on " << std::ctime(&time);
  std::string ts_str = timestamp.str();
  ts_str = std::strtok(&ts_str[0],"\n"); // Remove the \n appended by ctime

  const auto& filename = file_specs.filename;

  set_attribute<std::string>(filename,"case",p.get<std::string>("caseid","NONE"));
  set_attribute<std::string>(filename,"source","E3SM Atmosphere Model (EAMxx)");
  set_attribute<std::string>(filename,"eamxx_version",EAMXX_VERSION);
  set_attribute<std::string>(filename,"git_version",p.get<std::string>("git_version",EAMXX_GIT_VERSION));
  set_attribute<std::string>(filename,"hostname",p.get<std::string>("hostname","UNKNOWN"));
  set_attribute<std::string>(filename,"username",p.get<std::string>("username","UNKNOWN"));
  set_attribute<std::string>(filename,"atm_initial_conditions_file",p.get<std::string>("initial_conditions_file","NONE"));
  set_attribute<std::string>(filename,"topography_file",p.get<std::string>("topography_file","NONE"));
  set_attribute<std::string>(filename,"contact","e3sm-data-support@llnl.gov");
  set_attribute<std::string>(filename,"institution_id","E3SM-Projet");
  set_attribute<std::string>(filename,"realm","atmos");
  set_attribute<std::string>(filename,"history",ts_str);
  set_attribute<std::string>(filename,"Conventions","CF-1.8");
  set_attribute<std::string>(filename,"product",e2str(file_specs.ftype));
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
  m_atm_logger->info("           Filename prefix: " + m_filename_prefix);
  m_atm_logger->info("                    Run t0: " + m_run_t0.to_string());
  m_atm_logger->info("                   Case t0: " + m_case_t0.to_string());
  m_atm_logger->info("              Reference t0: " + m_output_control.last_write_ts.to_string());
  m_atm_logger->info("         Is Restart File ?: " + bool_to_string(m_is_model_restart_output));
  m_atm_logger->info("        Is Restarted Run ?: " + bool_to_string(m_is_restarted_run));
  m_atm_logger->info("            Averaging Type: " + e2str(m_avg_type));
  m_atm_logger->info("          Output Frequency: " + std::to_string(m_output_control.frequency) + " " + m_output_control.frequency_units);
  switch (m_output_file_specs.storage.type) {
    case NumSnaps:
      m_atm_logger->info("             File Capacity: " + std::to_string(m_output_file_specs.storage.max_snapshots_in_file) + " snapshots");  // TODO: add "not set" if the value is -1
      break;
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
