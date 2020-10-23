#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include <iostream>
#include <fstream>

#include "scream_config.h"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include "share/io/scream_scorpio_interface.hpp"

#include "share/field/field_repository.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_identifier.hpp"

#include "share/grid/grids_manager.hpp"

/*  The AtmosphereOutput class handles an output stream in SCREAM.
 *  Typical usage is to register an AtmosphereOutput object with the OutputManager, see output_manager.hpp
 *
 *  Similar to the typical AtmosphereProcess, output streams have a init, run and finalize routines.
 *  These routines are called during the respective steps of the AD call, i.e. ad.init(), ad.run() and ad.finalize()/
 *
 *  Each AtmosphereOutput instance handles one output stream.
 *
 *  At construction time ALL output instances require
 *  1. an EKAT comm group and
 *  2. a EKAT parameter list
 *  3. a shared pointer to the field repository
 *  4. a shared pointer to the grids manager
 * 
 *  The EKAT parameter list must contain all of the top level output control information and follows the format:
 *  ------
 *  FILENAME: STRING
 *  AVERAGING TYPE: STRING
 *  GRID: STRING                   (optional)
 *  FREQUENCY:
 *    OUT_N: INT
 *    OUT_OPTION: STRING
 *    OUT_MAX_STEPS: INT
 *  FIELDS:
 *    Number of Fields: INT
 *    field 1: STRING
 *    ...
 *    field N: STRING
 *  restart_hist_N: INT            (optional)
 *  restart_hist_OPTION: STRING    (optional)
 *  RESTART FILE: BOOL             (optional)
 *  -----
 *  where,
 *  FILENAME is a string of the filename suffix.  TODO: change this to a casename associated with the whole run.
 *  AVERAGING TYPE is a string that describes which type of output, current options are:
 *    instant - no averaging, output each snap as is.
 *    average - average of the field over some interval described in frequency section.
 *    min     - minimum value of the field over some time interval
 *    max     - maximum value of the field over some time interval
 *  GRID is a string describing which grid to write on, currently the only option is the default Physics.
 *  FREQUENCY is a subsection that controls the frequency parameters of output.
 *    OUT_N is an integer of the frequency of output given the units from OUT_OPTION
 *    OUT_OPTION is a string for the units of output frequency, examples would be "Steps", "Months", "Years", "Hours", etc.
 *    OUT_MAX_STEPS is an integer of the maximum number of steps that can exist on a single file (controls files getting too big).
 *  FIELDS is a subsection that lists all the fields in this output stream.
 *    Number of Fields is an integer that specifies the number of fields in this output stream.
 *    field 1,...,field N is a list of each field in the output stream by name in field manager.
 *  restart_hist_N is an optional integer parameter that specifies the frequenct of restart history writes.
 *  restart_hist_OPTION is an optional string parameter for the units of restart history output.
 *  RESTART FILE is an optional boolean parameter that specifies if this output stream is a restart output, which is treated differently.
 *
 *  Usage of this class is to create an output file, write data to the file and close the file.
 *  This class keeps a running copy of data for all output fields locally to be used for the different averaging flags.
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */

/* There are a number of things that need to be cleaned up - in no particular order -
 *--1) Gather dimension data from grid manager or field repo.  We may need to use field repo since grid manager
 *     is unaware of say the number of components, for example.  --DONE--
 *--2) Clean-up init step, perhaps with subfunctions.  --DONE--
 *--3) Gather degree's of freedom from grid manager. Need to know the global-id index for
 *     each value in a view, in order if flattened.  --DONE--
 *  4) Handle MPI communication to allow for PIO_STRIDE (i.e. fewer pio ranks than total ranks).
 *  5) Make it so that all dimensions that are registered are automatically registered as variables and the values written (only once, not with time).  Probably need this information from the grid manager(?)
 *--6) Write Restart output.  Extra-savy would be to use DAG from AD to determine which fields were essential and just write those.  #3
 *-6a) Also need restart history, in case you stop mid-out frequency and need to pick it up in the next run.
 *--7) Create an AtmosphereInput class built on these same concepts.  Can we make a master SCORPIO class with output and input as sub-classes?  --DONE--
 *--8) Create average, min, max and instantaneous options for output.  --DONE--
 *  9) Assure that IO is compatible with PACKing of variables.
 * 10) Add control of netCDF header information to the YAML file and this class.
 * 11) Better naming convention for output files?  Currently just keeps a counter, but in EAM the date/time range is used.
 * 12) Better handling of output frequency.  Need to a) gather time info from driver and b) parse that info for frequency units.
 * 13) Expand compatability of IO class to handle Dynamics grids.  This will require a different approach to "set_dofs"
 *-14) Hook up the AD for the unit test.  --DONE--
 *-15) SCORPIO in stand-alone build.
 * 16) When dynamics is ready, write output from a dynamics run.
 * 17) Remove the F90 layer and call PIO C routines directly.
 * 18) Add option to restart to write a restart at finalization, maybe if OUT_N = -1
 *-19) take advantage of shared pointers field_repo and gridsmanager to streamline init, run and finalize interface.
 * 20) Create a single IO Master class with AtmosphereOutput and AtmosphereInput as sub-classes.
 * 21) Output averaging types could be individual instances of the same class providing flexibility for future types to be incorporated.  See Jeff's work on field property tests for example.
 * 22) Currently the restart history frequency is set to the restart frequency.  Should maybe fix this so that if the restart history field is included in the parameter list for the output stream than that value is taken instead of the restart frequency values.
 * 23) Should restart output just be a sub-class of the AtmosphereOutput class?
 */

namespace scream
{

class AtmosphereOutput 
{
public:
  using value_type     = Real;          // Right not hard-coded to only Real type variables.  TODO: make this more general?
  using dofs_list_type = KokkosTypes<HostDevice>::view_1d<long>;
  using view_type      = typename KokkosTypes<HostDevice>::template view<value_type*>;
  using field_type     = Field<value_type>;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params, 
                   const std::shared_ptr<const FieldRepository<Real>>& repo,
                   const std::shared_ptr<const GridsManager>& gm)
  {
    m_comm       = comm;
    m_params     = params;
    m_field_repo = repo;
    m_gm         = gm;
  }

  // Main Functions
  void init();
  // Allow run to be called with a timestamp input (typical runs) and a floating point value (used for unit tests)
  void run(const Real time);
  void run(const util::TimeStamp& time);
  void finalize();

  // Helper Functions
  void check_status();
  std::map<std::string,Int> get_status() const { return m_status; }

protected:
  // Internal functions
  void add_identifier(const std::string name);
  void register_variables(const std::string filename);
  void set_degrees_of_freedom(const std::string filename);
  void register_views();
  void new_file(const std::string filename);
  void run_impl(const Real time, const std::string time_str);  // Actual run routine called by outward facing "run"
  // Internal variables
  ekat::ParameterList                          m_params;
  ekat::Comm                                   m_comm;
  std::shared_ptr<const FieldRepository<Real>> m_field_repo;
  std::shared_ptr<const GridsManager>          m_gm;
  
  // Main output control data
  std::string m_casename;
  std::string m_avg_type;
  std::string m_grid_name;
  // Frequency of output control
  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;
  // How individual columns are distributed across MPI Ranks
  Int m_total_dofs;
  Int m_local_dofs;
  // Restart history control
  Int m_restart_hist_n;
  std::string m_restart_hist_option;
  // Internal maps to the output fields, how the columns are distributed, the file dimensions and the global ids.
  std::map<std::string,field_type>       m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  dofs_list_type                         m_gids;
  // Local and global views of each field to be used for "averaging" output and writing to file.
  std::map<std::string,view_type>        m_view_global;
  std::map<std::string,view_type>        m_view_local;

  // Manage when files are open and closed, and what type of file I am writing.
  bool m_is_init = false;
  bool m_is_restart_hist = false;  //TODO:  If instead we rely on the timestamp to determine how many steps are represented in the averaging value, or maybe filename, we won't need this.
  bool m_is_restart = false;

  // Helper map to monitor an output stream's status.  Most used fields are Snaps and Avg Count which are
  // used to monitor if a new file is needed and if output should be written according to the frequency settings.
  std::map<std::string,Int> m_status = {
                                  {"Init",          0},  // Records the number of files this output stream has managed
                                  {"Run",           0},  // Total number of times "Run" has been called
                                  {"Finalize",      0},  // Total number of times "Finalize" has been called (should always be 1)
                                  {"Snaps",         0},  // Total number of timesnaps saved to the currently open file.
                                  {"Avg Count",     0},  // Total number of timesnaps that have gone by since the last time output was written.
                                       }; 
private:

}; // Class AtmosphereOutput

// ====================== IMPLEMENTATION ===================== //
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::init()
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
    // Find the <FieldIdentifier,Field> pair that corresponds with this variable name on the appropriate grid 
    add_identifier(var_name);
  }

  // Now that the fields have been gathered register the local and global views which will be used to determine output data to be written.
  register_views();

} // init
/* ---------------------------------------------------------- */
/* Overload the run routine to accept a TimeStamp or floating point for the input time */
inline void AtmosphereOutput::run(const util::TimeStamp& time)
{
  // In case it is needed for the output filename, parse the current timesnap into an appropriate string
  std::string time_str = time.to_string();
  std::replace( time_str.begin(), time_str.end(), ' ', '.');
  time_str.erase(std::remove( time_str.begin(), time_str.end(), ':'), time_str.end());
  // Pass the time in seconds and as a string to the run routine.
  run_impl(time.get_seconds(),time_str);
}
/*-----*/
inline void AtmosphereOutput::run(const Real time)
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
inline void AtmosphereOutput::run_impl(const Real time, const std::string time_str) 
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
  std::string filename = m_casename+"."+m_avg_type+"."+m_out_units+"_x"+std::to_string(m_out_frequency)+"."+time_str;
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
  //
  if (is_write) 
  {
    // If we closed the file in the last write because we reached max steps, or this is a restart history file,
    // we need to create a new file for writing.
    if( !is_typical or !m_is_init ) { new_file(filename) };
    if( !m_is_init and is_typical ) { m_is_init=true; }

    pio_update_time(filename,time); // Universal scorpio command to set the timelevel for this snap.
    if (is_typical) { m_status["Snaps"] += 1; }  // Update the snap tally, used to determine if a new file is needed and only needed for typical output.
  }

  // Take care of updating and possibly writing fields.
  for (auto const& fmap : m_fields)
  {
    // Get all the info for this field.
    auto name   = fmap.first;
    auto g_view = m_view_global.at(name); 
    Int  f_len  = fmap.second.get_header().get_identifier().get_layout().size();
    view_type l_view;
    // The next two operations do not need to happen if the frequency of output is instantaneous.
    if (m_avg_type == "Instant")
    {
      l_view = m_view_global.at(name);
    }
    else // output type uses multiple steps.
    {
      l_view = m_view_local.at(name);
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
      grid_write_data_array(filename,name,m_dofs.at(name),l_view.data());
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
inline void AtmosphereOutput::
finalize() 
{
  using namespace scream;
  using namespace scream::scorpio;

  // Nothing to do at the moment, but keep just in case future development needs a finalization step

  m_status["Finalize"] += 1;
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::check_status()
{
  printf("IO Status for Rank %5d, File - %.40s: (Init: %2d), (Run: %5d), (Finalize: %2d)\n",
                   m_comm.rank(),m_casename.c_str(),m_status["Init"],m_status["Run"],m_status["Finalize"]);
} // check_status
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::add_identifier(const std::string name)
{
/* 
 * add_identifier routine adds a new field identifier to the local map of fields that will be used by this class for IO.
 * INPUT:
 *   field_repo: is a pointer to the field_repository for this simulation.
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  auto f = m_field_repo->get_field(name, m_grid_name);
  auto fid = f.get_header().get_identifier();
  m_fields.emplace(name,f);
  // check to see if all the dims for this field are already set to be registered.
  for (int ii=0; ii<fid.get_layout().rank(); ++ii)
  {
    // check tag against m_dims map.  If not in there, then add it.
    auto& tag = fid.get_layout().tags()[ii];
    auto tag_loc = m_dims.find(tag2string(tag));
    if (tag_loc == m_dims.end()) 
    { 
      m_dims.emplace(std::make_pair(tag2string(tag),fid.get_layout().dim(ii)));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag2string(tag))==fid.get_layout().dim(ii),
        "Error! Dimension " + tag2string(tag) + " on field " + name + " has conflicting lengths");
    }
  }
} // add_identifier
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::register_views()
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto fmap : m_fields)
  {
    auto  name = fmap.first;
    // If the "averaging type" is instant then just need a ptr to the view.
    auto view_d = fmap.second.get_view();
    auto view_h = Kokkos::create_mirror_view( view_d );
    m_view_global.emplace(name,view_h);
    if (m_avg_type != "Instant") {
    // If anything other than instant we need a local copy that can be updated.
      view_type view_copy_d("",view_d.extent(0));
      auto view_copy = Kokkos::create_mirror_view( view_copy_d );
      m_view_local.emplace(name,view_copy);
    }
  }
}
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::register_variables(const std::string filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto fmap : m_fields)
  {
    auto  name = fmap.first;
    auto& fid  = fmap.second.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    for (auto& tag_ii : fid.get_layout().tags()) {
      io_decomp_tag += "-" + tag2string(tag_ii); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag2string(tag_ii)); // Add dimensions string to vector of dims.
    }
    io_decomp_tag += "-time";  // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger?  Should we register dimension variables (such as ncol and lat/lon) elsewhere in the dimension registration?  These won't have time.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.
    register_variable(filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
  }
  // Finish by registering time as a variable.  TODO: Should this really be something registered during the reg. dimensions step? 
  register_variable(filename,"time","time",1,{"time"},  PIO_REAL,"time");
  if (m_is_restart_hist) { register_variable(filename,"avg_count","avg_count",1,{"cnt"}, PIO_INT, "cnt"); }
} // register_variables
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::set_degrees_of_freedom(const std::string filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and set dof.
  for (auto fmap : m_fields)
  {
    auto  name = fmap.first;
    auto& fid  = fmap.second.get_header().get_identifier();
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
inline void AtmosphereOutput::new_file(const std::string filename)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Register new netCDF file for output.
  register_outfile(filename);

  // Register dimensions with netCDF file.
  for (auto it : m_dims)
  {
    if(it.first == "COL") { it.second=m_total_dofs;}  // Note: This is because only cols are decomposed over mpi ranks.  In this case still need max number of cols.
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
} //namespace scream
#endif // SCREAM_SCORPIO_OUTPUT_HPP
