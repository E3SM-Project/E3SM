#include "share/io/scorpio_output.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "share/io/scorpio_input.hpp"

#include "ekat/util/ekat_string_utils.hpp"

#include <numeric>
#include <fstream>

namespace scream
{

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const FieldManager<Real>>& field_mgr,
                  const bool is_restarted_run, const bool is_model_restart_output)
 : m_comm      (comm)
 , m_params    (params)
 , m_field_mgr (field_mgr)
 , m_is_model_restart_output (is_model_restart_output)
 , m_is_restarted_run (is_restarted_run)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (m_field_mgr, "Error! Invalid field manager pointer.\n");

  m_grid = m_field_mgr->get_grid();
  EKAT_REQUIRE_MSG(m_grid->name()=="Physics" || m_grid->name()=="Physics GLL",
      "Error! I/O only supports output on a Physics or Physics GLL grid.\n");

  EKAT_REQUIRE_MSG(m_comm.size()<=m_grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // Setup I/O structures
  init ();
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::init()
{
  // Parse the parameters that controls this output instance.
  m_casename = m_params.get<std::string>("Casename");

  m_max_snapshots_per_file  = m_params.get<int>("Max Snapshots Per File");
  m_out_frequency           = m_params.sublist("Output").get<int>("Frequency");
  m_out_frequency_units     = m_params.sublist("Output").get<std::string>("Frequency Units");

  auto avg_type = m_params.get<std::string>("Averaging Type");
  m_avg_type = ekat::upper_case(avg_type);
  auto valid_avg_types = {"INSTANT", "MAX", "MIN", "AVERAGE"};
  EKAT_REQUIRE_MSG (ekat::contains(valid_avg_types,m_avg_type),
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

  m_has_restart_data = (m_avg_type!="INSTANT");

  if (m_has_restart_data) {
    // This avg_type needs to save  some info in order to restart the output.
    // E.g., we might save 30-day avg value for field F, but due to job size
    // break the run into three 10-day runs. We then need to save the state of
    // our averaging in a "restart" file (e.g., the current avg and avg_count).
    // If checkpoint_freq=0, we do not perform checkpointing of the output.
    // NOTE: we may not need to *read* restart data (if this is not a restart run),
    //       but we might still want to save it. Viceversa, we might not want
    //       to save it, but still want to read it (a restarted run, where we
    //       don't want to save additional restart files).
    m_is_restarted_run     &= m_params.sublist("Restart").get("Perform Restart",true);
    m_checkpoint_freq       = m_params.sublist("Checkpointing").get("Frequency",0);
    m_checkpoint_freq_units = m_params.sublist("Checkpointing").get("Frequency Units",m_out_frequency_units);
  }

  // For each output field, ensure its dimensions are registered in the pio file
  m_fields = m_params.get<std::vector<std::string>>("Fields");
  for (const auto& var_name : m_fields) {
    register_dimensions(var_name);
  }

  // Now that the fields have been gathered register the local views which will be used to determine output data to be written.
  register_views();

  // If this is "normal" output of a restarted run, and avg_type requires restart data,
  // then we have to open the history restart file and read its data.
  // The user can skip this (which means the Output averaging would start from scratch)
  // by setting "Restart History" in the "Output Restart" to false.
  if (m_is_restarted_run && m_has_restart_data) {

    // TODO: add bool support to ekat::Comm MPI types
    bool found = false;
    std::string filename,rpointer_content;
    if (m_comm.am_i_root()) {
      std::ifstream rpointer_file;
      rpointer_file.open("rpointer.atm");
      auto hist_restart_testname = m_params.sublist("Restart").get("Casename",m_casename);
      auto testname = compute_filename_root(hist_restart_testname);
      while (rpointer_file >> filename) {
        rpointer_content += filename + "\n";
        if (filename.find(testname) != std::string::npos) {
          found = true;
          break;
        }
      }
    }

    // Sanity check
    m_comm.broadcast(&found,1,0);

    // If the history restart file is not found, we must HOPE that it is because
    // the last model restart step coincided with a model output step, in which case
    // a restart history file is not written.
    // TODO We have NO WAY of checking this from this class, but the OutputManager
    //      *might* be able to figure out if this is the case, and perhaps set
    //      "Restart History" to false in the input parameter list.
    //      For now, simply print a warning.
    if (found) {
      // Have the root rank communicate the nc filename
      int filename_size = filename.size();
      m_comm.broadcast(&filename_size,1,m_comm.root_rank());
      filename.resize(filename_size);
      m_comm.broadcast(&filename.front(),filename_size,m_comm.root_rank());

      // Create an input stream on the fly, and init averaging data
      ekat::ParameterList res_params("Input Parameters");
      res_params.set<std::string>("Filename",filename);
      res_params.set("Fields",m_fields);

      AtmosphereInput hist_restart (m_comm,res_params,m_grid,m_host_views_1d,m_layouts);
      hist_restart.read_variables();
      m_nsteps_since_last_output = hist_restart.read_int_scalar("avg_count");
      hist_restart.finalize();
    } else {
      if (m_comm.am_i_root()) {
        auto hist_restart_testname = m_params.sublist("Restart").get("Casename",m_casename);
        auto testname = compute_filename_root(hist_restart_testname);
        printf ("WARNING! No restart file found in the rpointer for case\n"
                "        %s\n"
                "   We *assume* this is because the last model restart write step\n"
                "   coincided with a model output step, so no history restart file\n"
                "   was needed. So far, we cannot check this, so we simply cross our fingers...\n",
                testname.c_str());
        printf ("rpointer.atm content:%s\n",rpointer_content.c_str());
      }
    }
  }

} // init
/* ---------------------------------------------------------- */
/* Overload the run routine to accept a TimeStamp or floating point for the input time */
void AtmosphereOutput::run(const util::TimeStamp& time)
{
  // In case it is needed for the output filename, parse the current timesnap into an appropriate string
  std::string time_str = time.to_string();
  std::replace( time_str.begin(), time_str.end(), ' ', '.');
  time_str.erase(std::remove( time_str.begin(), time_str.end(), ':'), time_str.end());
  // Pass the time in seconds and as a string to the run routine.
  run_impl(time.get_seconds(),time_str);
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::finalize() 
{
  // Safety check: we should not quit the simulation with an open file.
  // TODO: this might need to be adjusted once we test 2+ snapshots per file
  EKAT_REQUIRE_MSG (not m_is_output_file_open,
      "Error! AtmosphereOutput::finalize() was called while the output file was still open.\n");
} // finalize
/*-----*/
void AtmosphereOutput::run_impl(const Real time, const std::string& time_str) 
{
  using namespace scream::scorpio;

  ++m_nsteps_since_last_output;
  ++m_nsteps_since_last_checkpoint;
  // For the RUN step we always update the local views of each field to reflect the most recent step.
  // Following the update we have two courses of action:
  // 1. Do nothing else, this means that the frequency of output doesn't correspond with this step.
  // 2. Write output.
  //   a. This is either "normal" output or "restart" output.
  //   b. In the case of normal output we also reset the average counter.
  //   c. A model restart is just a normal output (with avg=Instant), so that fits under this category.
  // The other kind of output is the history restart.  This allows a restart run to also have
  // a consistent set of history. E.g., with avg_type=Average, and the avg window was not yet completed
  // when the model restart was written, we need to know what was the current value of avg (and avg_count)
  // at the restart point.
  // 1. A restart history is not necessary for,
  //   a. Instantaneous output streams.
  //   b. When also the model output was written during the step where model restart was written.
  //      In other words, when the average counter is 0.
  // Final point: a model output file, a model restart file, and an history restart file can all
  // be distinguished by the suffix of the generated file:
  //   '.nc': model output
  //   '.r.nc': model restart
  //   '.rhist.nc': history restart

  // Check to see if output is expected and what kind.
  // NOTE: for checkpoint, don't do m_nsteps_since_last_checkpoint==m_checkpoint_freq, since
  //   - you may have skipped a checkpoint since it coincided with output step
  //   - you can't zero m_nsteps_since_last_checkpoint if you do output, since that might
  //     get you out of sync with the model restart output (rhist should go hand-in-hand with
  //     model restart). E.g., if model restart (and history restart) every 3 days, but the
  //     model output is every 7 days, over a 30-day time horizon you would have:
  //       - model restarts at days 3,6,9,12,15,18,21,24,27,30.
  //       - model outputs at days 7,14,21,28.
  //       - history restarts at days 3,6,9,12,15,18,24,27,30, but not 21, cause model
  //         restart coincided with model output.
  //     If you zero nsteps_since_last_checkout at the 1st output, you would get history
  //     restarts all messed up, at day 10,13, ...
  const bool is_output_step = (m_nsteps_since_last_output == m_out_frequency);
  const bool is_checkpoint_step = !m_is_model_restart_output && m_checkpoint_freq>0 &&
                                   (m_nsteps_since_last_checkpoint % m_checkpoint_freq)==0 &&
                                  !is_output_step;

  // Output or checkpoint are both steps where we need to call scorpio
  const bool is_write_step = is_output_step || is_checkpoint_step;
  std::string filename;
  if (is_write_step) {
    filename = compute_filename_root(m_casename) + "." + time_str;
    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (m_is_model_restart_output) {
      filename+=".r.nc";
      if (m_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::out | std::ofstream::trunc);  // Open rpointer file and clear contents
        rpointer << filename << std::endl;
      }
    } else if (is_checkpoint_step) {
      filename+=".rhist.nc";
      if (m_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append the restart hist file information
        rpointer << filename << std::endl;
      }
    } else {
      filename += ".nc";
    }

    // If it's a checkpoint file, or if there's no file open for the output, we need to open the pio file.
    if( is_checkpoint_step or !m_is_output_file_open) {
      new_file(filename);
      if (is_checkpoint_step) { 
        set_int_attribute_c2f (filename.c_str(),"avg_count",m_nsteps_since_last_output);
      } else {
        m_is_output_file_open = true;
      }
    }

    // Set the timelevel for this snap in the pio file.
    pio_update_time(filename,time);
    if (is_output_step) {
      // We're adding one snapshot to the file
      ++m_num_snapshots_in_file;
    }
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields) {
    // Get all the info for this field.
    const auto  field = m_field_mgr->get_field(name);
    const auto& layout = m_layouts.at(name);
    const auto& dims = layout.dims();
    const auto  rank = layout.rank();

    // Safety check: make sure that the field was written at least once before using it.
    EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Output field '" + name + "' has not been initialized yet\n.");

    // Make sure host data is up to date
    field.sync_to_host();

    // Manually update the 'running-tally' views with data from the field,
    // by combining new data with current avg values.
    // NOTE: the running-tally is not a tally for Instant avg_type.
    auto data = m_host_views_1d.at(name).data();
    switch (rank) {
      case 1:
      {
        auto new_view_1d = field.get_view<const Real*,Host>();
        auto avg_view_1d = view_Nd_host<1>(data,dims[0]);

        for (int i=0; i<dims[0]; ++i) {
          combine(new_view_1d(i), avg_view_1d(i));
        }
        break;
      }
      case 2:
      {
        auto new_view_2d = field.get_view<const Real**,Host>();
        auto avg_view_2d = view_Nd_host<2>(data,dims[0],dims[1]);
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            combine(new_view_2d(i,j), avg_view_2d(i,j));
        }}
        break;
      }
      case 3:
      {
        auto new_view_3d = field.get_view<const Real***,Host>();
        auto avg_view_3d = view_Nd_host<3>(data,dims[0],dims[1],dims[2]);
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              combine(new_view_3d(i,j,k), avg_view_3d(i,j,k));
        }}}
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
    }

    if (is_write_step) {
      grid_write_data_array(filename,name,data);
    }
  }

  // Finish up any updates to output file and snap counter.
  if (is_write_step) {
    sync_outfile(filename);
    // If snaps equals max per file, close this file and set flag to open a new one next write step.
    if (is_output_step) {
      if (m_num_snapshots_in_file == m_max_snapshots_per_file) {
        m_num_snapshots_in_file = 0;

        // This file is "full". Close it.
        eam_pio_closefile(filename); 

        // Make sure a new output file will be created at the next output step.
        m_is_output_file_open = false;
      }
      // Zero out the Avg Count count now that snap has been written.
      m_nsteps_since_last_output = 0;

      // Since we saved, 
    } else {
      // A checkpoint step, close the file.
      eam_pio_closefile(filename); 

      // Zero out the checkpoint step counter
      m_nsteps_since_last_checkpoint = 0;
    }
  }
} // run_impl

/* ---------------------------------------------------------- */
void AtmosphereOutput::register_dimensions(const std::string& name)
{
/* 
 * Checks that the dimensions associated with a specific variable will be registered with IO file.
 * INPUT:
 *   field_manager: is a pointer to the field_manager for this simulation.
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  using namespace scorpio;

  // Store the field layout
  const auto& fid = m_field_mgr->get_field(name).get_header().get_identifier();
  const auto& layout = fid.get_layout();
  m_layouts.emplace(name,layout);

  // Now check taht all the dims of this field are already set to be registered.
  for (int i=0; i<layout.rank(); ++i) {
    // check tag against m_dims map.  If not in there, then add it.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    const auto tag_name = get_nc_tag_name(tags[i],dims[i]);
    auto tag_loc = m_dims.find(tag_name);
    if (tag_loc == m_dims.end()) {
      int tag_len = 0;
      if(e2str(tags[i]) == "COL") {
        // Note: This is because only cols are decomposed over mpi ranks. 
        //       In this case, we need the GLOBAL number of cols.
        tag_len = m_grid->get_num_global_dofs();
      } else {
        tag_len = layout.dim(i);
      }
      m_dims.emplace(std::make_pair(get_nc_tag_name(tags[i],dims[i]),tag_len));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag_name)==dims[i] or e2str(tags[i])=="COL",
        "Error! Dimension " + tag_name + " on field " + name + " has conflicting lengths");
    }
  }
} // register_dimensions
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_views()
{
  // Cycle through all fields and register.
  for (auto const& name : m_fields) {
    auto field = m_field_mgr->get_field(name);

    // These local views are really only needed if the averaging time is not 'Instant',
    // to store running tallies for the average operation. However, we create them
    // also for Instant avg_type, for simplicity later on.

    // If we have an 'Instant' avg type, we can alias the tmp views with the
    // host view of the field, provided that the field does not have padding,
    // and that it is not a subfield of another field (or else the view
    // would be strided).
    bool can_alias_field_view =
        m_avg_type=="Instant" &&
        field.get_header().get_alloc_properties().get_padding()==1 &&
        field.get_header().get_parent().expired();

    const auto size = m_layouts.at(name).size();
    if (can_alias_field_view) {
      // Alias field's data, to save storage.
      m_host_views_1d.emplace(name,view_1d_host(field.get_internal_view_data<Host>(),size));
    } else {
      // Create a local host view.
      m_host_views_1d.emplace(name,view_1d_host("",size));
    }
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_variables(const std::string& filename)
{
  using namespace scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    const auto& layout = fid.get_layout();
    for (int i=0; i<fid.get_layout().rank(); ++i) {
      const auto tag_name = get_nc_tag_name(layout.tag(i), layout.dim(i));
      io_decomp_tag += "-" + tag_name; // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }
    // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger? 
    // Should we register dimension variables (such as ncol and lat/lon) elsewhere
    // in the dimension registration?  These won't have time. 
    io_decomp_tag += "-time";
    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to fully C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.

     // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.
    register_variable(filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);
  }
  // Finish by registering time as a variable.
  // TODO: Should this really be something registered during the reg. dimensions step? 
  register_variable(filename,"time","time",1,{"time"},  PIO_REAL,"time");
} // register_variables
/* ---------------------------------------------------------- */
std::vector<int> AtmosphereOutput::get_var_dof_offsets(const FieldLayout& layout)
{
  std::vector<int> var_dof(layout.size());

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // These are not the dofs global ids (which are just labels, and can be whatever,
  // and in fact are not even contiguous when Homme generates the dof gids).
  // So, if the returned vector is {2,3,4,5}, it means that the 4 dofs on this rank
  // correspond to the 3rd,4th,5th, and 6th dofs globally.
  if (layout.has_tag(ShortFieldTagsNames::COL)) {
    const int num_cols = m_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smalle.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    int col_size = layout.size() / num_cols;

    // Compute the number of columns owned by all previous ranks.
    int offset = 0;
    m_comm.scan(&num_cols,&offset,1,MPI_SUM);
    offset -= num_cols;

    // Compute offsets of all my dofs
    std::iota(var_dof.begin(), var_dof.end(), offset*col_size);
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_degrees_of_freedom(const std::string& filename)
{
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // Cycle through all fields and set dof.
  for (auto const& name : m_fields) {
    auto field = m_field_mgr->get_field(name);
    const auto& fid  = field.get_header().get_identifier();
    auto var_dof = get_var_dof_offsets(fid.get_layout());
    set_dof(filename,name,var_dof.size(),var_dof.data());
    m_dofs.emplace(std::make_pair(name,var_dof.size()));
  }
  // Set degree of freedom for "time"
  int time_dof[1] = {0};
  set_dof(filename,"time",0,time_dof);

  /* TODO: 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
void AtmosphereOutput::new_file(const std::string& filename)
{
  using namespace scream::scorpio;

  // Register new netCDF file for output.
  register_outfile(filename);

  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    register_dimension(filename,it.first,it.first,it.second);
  }
  // Note: time has an unknown length. Setting its "length" to 0 tells the scorpio to
  // set this dimension as having an 'unlimited' length, thus allowing us to write
  // as many timesnaps to file as we desire.
  register_dimension(filename,"time","time",0);

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename);

  // Set the offsets of the local dofs in the global vector.
  set_degrees_of_freedom(filename);

  // Finish the definition phase for this file.
  eam_pio_enddef (filename); 
}

std::string AtmosphereOutput::compute_filename_root (const std::string& casename) const
{
  return casename + "." +
         m_avg_type + "." +
         m_out_frequency_units + "_x" +
         std::to_string(m_out_frequency);
}

} // namespace scream
