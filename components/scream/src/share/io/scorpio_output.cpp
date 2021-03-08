#include "share/io/scorpio_output.hpp"

namespace scream
{

// ====================== IMPLEMENTATION ===================== //
/* ---------------------------------------------------------- */
void AtmosphereOutput::init()
{
  using namespace scream;
  using namespace scream::scorpio;

  // Parse the parameters that controls this output instance.
  // See the comments at the top for more details.
  m_casename        = m_params.get<std::string>("FILENAME");
  m_avg_type        = m_params.get<std::string>("AVERAGING TYPE");
  m_grid_name       = m_params.get<std::string>("GRID","Physics");  // optional, default to Physics.
  auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");
  m_restart_hist_n  = m_params.get<Int>("restart_hist_N",0);  // optional, default to 0 for no output
  m_restart_hist_option = m_params.get<std::string>("restart_hist_OPTION","NONE"); // optional, default to NONE
  m_is_restart = m_params.get<bool>("RESTART FILE",false);  // optional, default to false

  // Gather data from grid manager:  In particular the global ids for columns assigned to this MPI rank
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with output grid! scorpio_output.hpp class only supports output on a Physics grid for now.\n");
  auto gids_dev = m_gm->get_grid(m_grid_name)->get_dofs_gids();
  m_gids = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids,gids_dev); 
  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  // int total_dofs;
  m_local_dofs = m_gids.size();
  MPI_Allreduce(&m_local_dofs, &m_total_dofs, 1, MPI_INT, MPI_SUM, m_comm.mpi_comm());
  EKAT_REQUIRE_MSG(m_comm.size()<=m_total_dofs,"Error, PIO interface only allows for the IO comm group size to be less than or equal to the total # of columns in grid.  Consider decreasing size of IO comm group.\n");

  // Create map of fields in this output with the field_identifier in the field repository.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    // Determine the variable name 
    std::string var_name = var_params.get<std::string>(ekat::strint("field",var_i+1));
    m_fields.push_back(var_name);
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
      printf("Warning! No restart history file found in rpointer file for %s, using current values in field repo\n",m_casename.c_str());
    }
    // Register rhist file as input and copy data to local views
    ekat::ParameterList res_params("Input Parameters");
    res_params.set<std::string>("FILENAME",filename);
    res_params.set<std::string>("GRID","Physics");
    res_params.set<bool>("RHIST",true);
    auto& f_list = res_params.sublist("FIELDS");
    Int fcnt = 1;
    f_list.set<Int>("Number of Fields",m_fields.size());
    for (auto name : m_fields)
    {
      f_list.set<std::string>("field "+std::to_string(fcnt),name);
      fcnt+=1;
    }
    input_type rhist_in(m_comm,res_params,m_field_repo,m_gm);
    rhist_in.init();
    for (auto name : m_fields)
    {
      auto l_view = rhist_in.pull_input(name);
      m_view_local.at(name) = l_view;
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
  using namespace scream;
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
    if( !m_is_init and is_typical ) { m_is_init=true; }

    pio_update_time(filename,time); // Universal scorpio command to set the timelevel for this snap.
    if (is_typical) { m_status["Snaps"] += 1; }  // Update the snap tally, used to determine if a new file is needed and only needed for typical output.
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields)
  {
    // Get all the info for this field.
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto view_d = field.get_view();
    auto g_view = Kokkos::create_mirror_view( view_d );
    Kokkos::deep_copy(g_view, view_d);
    Int  f_len  = field.get_header().get_identifier().get_layout().size();
    auto l_view = m_view_local.at(name);
    // The next two operations do not need to happen if the frequency of output is instantaneous.
    if (m_avg_type == "Instant")
    {
      Kokkos::deep_copy(l_view, g_view);
    }
    else // output type uses multiple steps.
    {
      // Update local view given the averaging type.  TODO make this a switch statement?
      if (m_avg_type == "Average")
      {
        for (int ii=0; ii<f_len; ++ii) {l_view(ii) = (l_view(ii)*(m_status["Avg Count"]-1) + g_view(ii))/(m_status["Avg Count"]);}
      } else if (m_avg_type == "Max")
      {
        for (int ii=0; ii<f_len; ++ii) {l_view(ii) = std::max(l_view(ii),g_view(ii));}
      } else if (m_avg_type == "Min")
      {
        for (int ii=0; ii<f_len; ++ii) {l_view(ii) = std::min(l_view(ii),g_view(ii));}
      } else {
        EKAT_REQUIRE_MSG(true, "Error! IO Class, updating local views, averaging type of " + m_avg_type + " is not supported.");
      }
    } // m_avg_type != "Instant"
    if (is_write) {
      auto l_dims = field.get_header().get_identifier().get_layout().dims();
      Int padding = field.get_header().get_alloc_properties().get_padding();
      grid_write_data_array(filename,name,l_dims,m_dofs.at(name),padding,l_view.data());
      if (is_typical) { 
        for (int ii=0; ii<f_len; ++ii) { l_view(ii) = g_view(ii); }  // Reset local view after writing.  Only for typical output.
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
  using namespace scream;
  using namespace scream::scorpio;

  // Nothing to do at the moment, but keep just in case future development needs a finalization step

  m_status["Finalize"] += 1;
} // finalize
void AtmosphereOutput::register_dimensions(const std::string& name)
{
/* 
 * Checks that the dimensions associated with a specific variable will be registered with IO file.
 * INPUT:
 *   field_repo: is a pointer to the field_repository for this simulation.
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  auto fid = m_field_repo->get_field(name, m_grid_name).get_header().get_identifier();
  // check to see if all the dims for this field are already set to be registered.
  for (int ii=0; ii<fid.get_layout().rank(); ++ii)
  {
    // check tag against m_dims map.  If not in there, then add it.
    auto& tag = fid.get_layout().tags()[ii];
    auto tag_loc = m_dims.find(e2str(tag));
    if (tag_loc == m_dims.end()) 
    { 
      Int tag_len = 0;
      if(e2str(tag) == "COL")
      {
        tag_len = m_total_dofs;  // Note: This is because only cols are decomposed over mpi ranks.  In this case still need max number of cols.
      } else {
        tag_len = fid.get_layout().dim(ii);
      }
      m_dims.emplace(std::make_pair(e2str(tag),tag_len));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(e2str(tag))==fid.get_layout().dim(ii) or e2str(tag)=="COL",
        "Error! Dimension " + e2str(tag) + " on field " + name + " has conflicting lengths");
    }
  }
} // register_dimensions
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_views()
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields)
  {
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto ts    = field.get_header().get_tracking().get_time_stamp();
    EKAT_REQUIRE_MSG (ts.is_valid(), "Error in output, field " + name + " has not been set in the field manager yet");
    // If the "averaging type" is instant then just need a ptr to the view.
    EKAT_REQUIRE_MSG (field.get_header().get_parent().expired(), "Error! Cannot deal with subfield, for now.");
    auto view_d = field.get_view();
    // Create a local copy of view to be stored by output stream.
    view_type_host view_copy("",view_d.extent(0));
    Kokkos::deep_copy(view_copy, view_d);
    m_view_local.emplace(name,view_copy);
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_variables(const std::string& filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields)
  {
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = field.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    for (auto& tag_ii : fid.get_layout().tags()) {
      io_decomp_tag += "-" + e2str(tag_ii); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(e2str(tag_ii)); // Add dimensions string to vector of dims.
    }
    io_decomp_tag += "-time";  // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger?  Should we register dimension variables (such as ncol and lat/lon) elsewhere in the dimension registration?  These won't have time.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.
    register_variable(filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
  }
  // Finish by registering time as a variable.  TODO: Should this really be something registered during the reg. dimensions step? 
  register_variable(filename,"time","time",1,{"time"},  PIO_REAL,"time");
  if (m_is_restart_hist) { register_variable(filename,"avg_count","avg_count",1,{"cnt"}, PIO_REAL, "cnt"); }
} // register_variables
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_degrees_of_freedom(const std::string& filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and set dof.
  for (auto const& name : m_fields)
  {
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = field.get_header().get_identifier();
    // bool has_cols = true;
    Int dof_len, n_dim_len, num_cols;
    // Total number of values represented by this rank for this field is given by the field_layout size.
    dof_len = fid.get_layout().size();
    // For a SCREAM Physics grid, only the total number of columns is decomposed over MPI ranks.  The global id (gid)
    // is stored here as m_gids.  Thus, for this field, the total number of dof's in the other dimensions (i.e. levels)
    // can be found by taking the quotient of dof_len and the length of m_gids.
    if (fid.get_layout().has_tag(FieldTag::Column))
    {
      num_cols = m_gids.size();
    } else {
      // This field is not defined over columns
      // TODO, when we allow for dynamics mesh this check will need to be adjusted for the element tag as well.
      num_cols = 1;
      // has_cols = false;
    }
    n_dim_len = dof_len/num_cols;
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global output indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    std::vector<Int> var_dof(dof_len);
    Int dof_it = 0;
    for (int ii=0;ii<num_cols;++ii)
    {
      for (int jj=0;jj<n_dim_len;++jj)
      {
        var_dof[dof_it] =  m_gids(ii)*n_dim_len + jj;
        ++dof_it;
      }
    }
    set_dof(filename,name,dof_len,var_dof.data());
    m_dofs.emplace(std::make_pair(name,dof_len));
  }
  // Set degree of freedom for "time"
  set_dof(filename,"time",0,0);
  if (m_is_restart_hist) { 
    Int var_dof[1] = {0};
    set_dof(filename,"avg_count",1,var_dof); 
   }
  /* TODO: 
   * Adjust DOF to accomodate packing for fields 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
void AtmosphereOutput::new_file(const std::string& filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Register new netCDF file for output.
  register_outfile(filename);

  // Register dimensions with netCDF file.
  for (auto it : m_dims)
  {
    register_dimension(filename,it.first,it.first,it.second);
  }
  register_dimension(filename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  if (m_is_restart_hist) { register_dimension(filename,"cnt","cnt",1); }
  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename);
  set_degrees_of_freedom(filename);

  // Finish the definition phase for this file.
  eam_pio_enddef  (filename); 

}
/* ---------------------------------------------------------- */
} // namespace scream
