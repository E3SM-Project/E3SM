#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"

#include <numeric>
#include "ekat/util/ekat_string_utils.hpp"

namespace scream
{

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const FieldManager<Real>>& field_mgr,
                  const std::shared_ptr<const GridsManager>& grid_mgr,
                  const bool read_restart_hist)
 : m_comm      (comm)
 , m_params    (params)
 , m_field_mgr (field_mgr)
 , m_grid_mgr  (grid_mgr)
 , m_read_restart_hist (read_restart_hist)
{
  EKAT_REQUIRE_MSG (m_field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (m_grid_mgr,  "Error! Invalid grid manager pointer.\n");
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::init()
{
  // Parse the parameters that controls this output instance.
  // See the comments at the top for more details.
  m_casename        = m_params.get<std::string>("FILENAME");
  m_grid_name       = m_params.get<std::string>("GRID","Physics");  // optional, default to Physics.

  const auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");

  m_restart_hist_n  = m_params.get<Int>("restart_hist_N",0);  // optional, default to 0 for no output
  m_restart_hist_option = m_params.get<std::string>("restart_hist_OPTION","NONE"); // optional, default to NONE
  m_is_restart = m_params.get<bool>("RESTART FILE",false);  // optional, default to false

  m_avg_type = m_params.get<std::string>("AVERAGING TYPE");
  auto avg_type_ci = ekat::upper_case(m_avg_type);

  EKAT_REQUIRE_MSG (
      avg_type_ci=="INSTANT" || avg_type_ci=="AVERAGE" || avg_type_ci=="MAX" || avg_type_ci=="MIN",
      "Error! Unsupported averaging type (" + m_avg_type + ").\n"
      "       Possible choices: Instant, Average, Max, Min.\n");

  // Gather data from grid manager:  In particular the global ids for columns assigned to this MPI rank
  EKAT_REQUIRE_MSG(m_grid_name=="Physics" || m_grid_name=="Physics GLL",
      "Error! I/O only supports output on a Physics or Physics GLL grid.\n");

  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  // int total_dofs;
  m_total_dofs = m_grid_mgr->get_grid(m_grid_name)->get_num_global_dofs();
  EKAT_REQUIRE_MSG(m_comm.size()<=m_total_dofs,
      "Error! PIO interface requires the size of the IO MPI group to be no greater\n"
      "       than the global number of columns. Consider decreasing the size of IO MPI group.\n");

  // Create map of fields in this output with the field_identifier in the field manager.
  m_fields = m_params.get<std::vector<std::string>>("FIELDS");
  for (const auto& var_name : m_fields) {
    /* Check that all dimensions for this variable are set to be registered */
    register_dimensions(var_name);
  }

  // Now that the fields have been gathered register the local views which will be used to determine output data to be written.
  register_views();

  // If this is a restart run that requires a restart history file read input here:
  if (m_read_restart_hist)
  {
    std::ifstream rpointer_file;
    rpointer_file.open("rpointer.atm");
    std::string filename;
    bool found = false;
    std::string testname = m_casename+"."+m_avg_type+"."+m_out_units+"_x"+std::to_string(m_out_frequency);
    while (rpointer_file >> filename)
    {
      if (filename.find(testname) != std::string::npos)
      {
        found = true;
        break;
      }
    }
    if ( not found )
    {
      printf("Warning! No restart history file found in rpointer file for %s, using current values in field manager\n",m_casename.c_str());
    }
    // Register rhist file as input and copy data to local views
    ekat::ParameterList res_params("Input Parameters");
    res_params.set<std::string>("FILENAME",filename);
    res_params.set<std::string>("GRID","Physics");
    res_params.set<bool>("RHIST",true);
    res_params.set("FIELDS",m_fields);

    using input_type     = AtmosphereInput;
    input_type rhist_in(m_comm,res_params,m_field_mgr,m_grid_mgr);
    rhist_in.init();
    for (auto name : m_fields)
    {
      auto l_view = rhist_in.pull_input(name);
      m_host_views_1d.at(name) = l_view;
    }
    auto avg_count = rhist_in.pull_input("avg_count");
    m_status["Avg Count"] = avg_count(0);
    rhist_in.finalize();
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
/*-----*/
void AtmosphereOutput::run(const Real time)
{
  // In case it is needed for the output filename, parse the current timesnap into an appropriate string
  // Convert time in seconds to a DD-HHMMSS string:
  const int ss = static_cast<int>(time);
  const int h =  ss / 3600;
  const int m = (ss % 3600) / 60;
  const int s = (ss % 3600) % 60;
  const std::string zero = "00";
  std::string time_str = (h==0 ? zero : std::to_string(h)) + (m==0 ? zero : std::to_string(m)) + (s==0 ? zero : std::to_string(s));
  // Pass the time in seconds and as a string to the run routine.
  run_impl(time,time_str);
}
/*-----*/
void AtmosphereOutput::run_impl(const Real time, const std::string& time_str) 
{
  using namespace scream::scorpio;

  m_status["Run"] += 1;
  m_status["Avg Count"] += 1;
  // For the RUN step we always update the local views of each field to reflect the most recent step.
  // Following the update we have two courses of action:
  // 1. Do nothing else, this means that the frequency of output doesn't correspond with this step.
  // 2. Write output.
  //   a. This is either typical output or restart output.
  //   b. In the case of typical output we also reset the average counter and the local view.
  //   c. A standard restart is just an instance of typical output, so that fits under this category.
  // The other kind of output is the restart history output.  This allows a restart run to also have
  // a consistent set if history outputs.
  // 1. A restart history is not necessary for,
  //   a. Instantaneous output streams.
  //   b. When the history has also be written this step.  In other words, when the average counter is 0.
  // Final point, typical output, a restart file and a restart history file can be distinguished by the
  // suffix, (.nc) is typical, (.r.nc) is a restart and (.rhist.nc) is a restart history file.

  // Check to see if output is expected and what kind.
  bool is_typical = (m_status["Avg Count"] == m_out_frequency);  // It is time to write output data.
  bool is_rhist   = (m_status["Avg Count"] == m_restart_hist_n) and !is_typical;  // It is time to write a restart history file.
  bool is_write = is_typical or is_rhist; // General flag for if output is written
  // Preamble to writing output this step
  std::string filename;
  if (is_write)
  {
    filename = m_casename+"."+m_avg_type+"."+m_out_units+"_x"+std::to_string(m_out_frequency)+"."+time_str;
    // Typical out can still be restart output if this output stream is for a restart file.  If it is a restart file it has a different suffix
    // and the filename needs to be added to the rpointer.atm file.
    if (m_is_restart) 
    { 
      filename+=".r";
      std::ofstream rpointer;
      rpointer.open("rpointer.atm",std::ofstream::out | std::ofstream::trunc);  // Open rpointer file and clear contents
      rpointer << filename + ".nc" << std::endl;
    }
    // If the output written will be to a restart history file than make sure the suffix is correct.
    if (is_rhist)
    { 
      filename+=".rhist"; 
      std::ofstream rpointer;
      rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append the restart hist file information
      rpointer << filename + ".nc" << std::endl;
      m_is_restart_hist = true;
    }
    filename += ".nc";
    // If we closed the file in the last write because we reached max steps, or this is a restart history file,
    // we need to create a new file for writing.
    if( !is_typical or !m_is_init ) 
    { 
      new_file(filename);
      if (is_rhist) 
      { 
        std::array<Real,1> avg_cnt = { (Real) m_status["Avg Count"] };
        grid_write_data_array(filename,"avg_count",avg_cnt.size(),avg_cnt.data());
      }
    }
    // Now the filename that is being stored in this object should be the appropriate file to be writing too.
    filename = m_filename;
    if( !m_is_init and is_typical ) { m_is_init=true; }

    pio_update_time(filename,time); // Universal scorpio command to set the timelevel for this snap.
    if (is_typical) { m_status["Snaps"] += 1; }  // Update the snap tally, used to determine if a new file is needed and only needed for typical output.
  }

  // Take care of updating and possibly writing fields.
  using view_2d_host = typename KokkosTypes<HostDevice>::view_2d<Real>;
  using view_3d_host = typename KokkosTypes<HostDevice>::view_3d<Real>;
  for (auto const& name : m_fields)
  {
    // Get all the info for this field.
    const auto  field = m_field_mgr->get_field(name);
    const auto& fap = field.get_header().get_alloc_properties();
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto rank = layout.rank();
    const auto num_reals = fap.get_alloc_size() / sizeof(Real);
    const auto last_dim = fap.get_last_extent();

    // Make sure host data is up to date
    field.sync_to_host();

    // This tmp view may contain garbage entries, if the field has padding
    view_1d_host hview_1d("",num_reals);
    switch (rank) {
      case 1:
      {
        auto view_1d = field.get_view<const Real*,Host>();
        Kokkos::deep_copy(hview_1d,view_1d);
        break;
      }
      case 2:
      {
        auto view_2d = field.get_view<const Real**,Host>();
        Kokkos::deep_copy(view_2d_host(hview_1d.data(),layout.dim(0),last_dim),view_2d);
        break;
      }
      case 3:
      {
        auto view_3d = field.get_view<const Real***,Host>();
        Kokkos::deep_copy(view_3d_host(hview_1d.data(),layout.dim(0),layout.dim(1),last_dim),view_3d);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
    }

    Int  f_len  = field.get_header().get_identifier().get_layout().size();
    auto l_view = m_host_views_1d.at(name);
    // It is not necessary to do any operations between local and global views if the frequency of output
    // is instantaneous, or if the Average Counter is 1 (meaning the beginning of a new record).
    // TODO: Question to address - This current approach will *not* include the initial conditions in the calculation of any of the
    // output metrics.  Do we want this to be the case? 
    if (m_avg_type == "Instant" || m_status["Avg Count"] == 1) {
      // Make sure that the global view is in fact valid before copying to local view.
      EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
          "Error in output, field " + name + " has not been initialized yet\n.");
      Kokkos::deep_copy(l_view, hview_1d);
    } else {
      // output type uses multiple snapshots.

      // Update local view given the averaging type.  TODO make this a switch statement?
      if (m_avg_type == "Average") {
        for (int ii=0; ii<f_len; ++ii) {
          l_view(ii) = (l_view(ii)*(m_status["Avg Count"]-1) + hview_1d(ii))/(m_status["Avg Count"]);
        }
      } else if (m_avg_type == "Max") {
        for (int ii=0; ii<f_len; ++ii) {
          l_view(ii) = std::max(l_view(ii),hview_1d(ii));
        }
      } else if (m_avg_type == "Min") {
        for (int ii=0; ii<f_len; ++ii) {
          l_view(ii) = std::min(l_view(ii),hview_1d(ii));
        }
      }
    } // m_avg_type != "Instant"

    if (is_write) {
      auto l_dims = field.get_header().get_identifier().get_layout().dims();
      Int padding = field.get_header().get_alloc_properties().get_padding();
      grid_write_data_array(filename,name,l_dims,m_dofs.at(name),padding,l_view.data());
      if (is_typical) { 
        for (int ii=0; ii<f_len; ++ii) {
          l_view(ii) = hview_1d(ii);
        }  // Reset local view after writing.  Only for typical output.
      }
    }
  }

  // Finish up any updates to output file and snap counter.
  if (is_write)
  {
    sync_outfile(filename);
    // If snaps equals max per file, close this file and set flag to open a new one next write step.
    if (is_typical)
    {
      if (m_status["Snaps"] == m_out_max_steps)
      {
        m_status["Snaps"] = 0;
        eam_pio_closefile(filename);
        m_is_init = false;
      }
      // Zero out the Avg Count count now that snap has been written.
      m_status["Avg Count"] = 0;
    }
    else  // must be that is_rhist=true
    { 
      eam_pio_closefile(filename); 
    }
  }
  // Reset flag for restart history write.
  m_is_restart_hist = false;

} // run
/* ---------------------------------------------------------- */
void AtmosphereOutput::finalize() 
{
  // Nothing to do at the moment, but keep just in case future development needs a finalization step

  m_status["Finalize"] += 1;
} // finalize
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

  auto fid = m_field_mgr->get_field(name).get_header().get_identifier();
  // check to see if all the dims for this field are already set to be registered.
  for (int ii=0; ii<fid.get_layout().rank(); ++ii) {
    // check tag against m_dims map.  If not in there, then add it.
    const auto& tags = fid.get_layout().tags();
    const auto& dims = fid.get_layout().dims();
    const auto tag_name = get_nc_tag_name(tags[ii],dims[ii]);
    auto tag_loc = m_dims.find(tag_name);
    if (tag_loc == m_dims.end()) {
      Int tag_len = 0;
      if(e2str(tags[ii]) == "COL") {
        tag_len = m_total_dofs;  // Note: This is because only cols are decomposed over mpi ranks.  In this case still need max number of cols.
      } else {
        tag_len = fid.get_layout().dim(ii);
      }
      m_dims.emplace(std::make_pair(get_nc_tag_name(tags[ii],dims[ii]),tag_len));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag_name)==fid.get_layout().dim(ii) or e2str(tags[ii])=="COL",
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

    // Create a local host view.
    const auto num_reals = field.get_header().get_alloc_properties().get_alloc_size() / sizeof(Real);
    m_host_views_1d.emplace(name,view_1d_host("",num_reals));
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
  if (m_is_restart_hist) {
    register_variable(filename,"avg_count","avg_count",1,{"cnt"}, PIO_REAL, "cnt");
  }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<Int> AtmosphereOutput::get_var_dof_offsets(const int dof_len, const bool has_cols)
{
  std::vector<Int> var_dof(dof_len);

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // These are not the dofs global ids (which are just labels, and can be whatever,
  // and in fact are not even contiguous when Homme generates the dof gids).
  // So, if the returned vector is {2,3,4,5}, it means that the 4 dofs on this rank
  // correspond to the 3rd,4th,5th, and 6th dofs globally.
  if (has_cols) {
    const int num_cols = m_grid_mgr->get_grid(m_grid_name)->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smalle.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    Int col_size = dof_len/num_cols;

    // Compute the number of columns owned by all previous ranks.
    Int offset = 0;
    m_comm.scan_sum(&num_cols,&offset,1);

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
    auto& fid  = field.get_header().get_identifier();
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global output indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    const bool has_col_tag = fid.get_layout().has_tag(COL);
    std::vector<Int> var_dof = get_var_dof_offsets(fid.get_layout().size(), has_col_tag);
    set_dof(filename,name,var_dof.size(),var_dof.data());
    m_dofs.emplace(std::make_pair(name,var_dof.size()));
  }
  // Set degree of freedom for "time"
  set_dof(filename,"time",0,0);
  if (m_is_restart_hist) { 
    Int var_dof[1] = {0};
    set_dof(filename,"avg_count",1,var_dof); 
   }
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
  m_filename = filename;

  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    register_dimension(filename,it.first,it.first,it.second);
  }
  register_dimension(filename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  if (m_is_restart_hist) { register_dimension(filename,"cnt","cnt",1); }
  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename);
  set_degrees_of_freedom(filename);

  // Finish the definition phase for this file.
  eam_pio_enddef (filename); 
}

} // namespace scream
