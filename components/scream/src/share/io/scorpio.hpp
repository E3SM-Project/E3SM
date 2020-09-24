#ifndef SCREAM_SCORPIO_HPP
#define SCREAM_SCORPIO_HPP

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
/* There are a number of things that need to be cleaned up - in no particular order -
 *--1) Gather dimension data from grid manager or field repo.  We may need to use field repo since grid manager
 *     is unaware of say the number of components, for example.  --DONE--
 *--2) Clean-up init step, perhaps with subfunctions.  --DONE--
 *--3) Gather degree's of freedom from grid manager. Need to know the global-id index for
 *     each value in a view, in order if flattened.  --DONE--
 *  4) Handle MPI communication to allow for PIO_STRIDE (i.e. fewer pio ranks than total ranks).
 *  5) Make it so that all dimensions that are registered are automatically registered as variables and the values written (only once, not with time).  Probably need this information from the grid manager(?)
 *  6) Write Restart output.  Extra-savy would be to use DAG from AD to determine which fields were essential and just write those.
 *  7) Create an AtmosphereInput class built on these same concepts.  Can we make a master SCORPIO class with output and input as sub-classes?
 *--8) Create average, min, max and instantaneous options for output.  --DONE--
 *  9) Assure that IO is compatible with PACKing of variables.
 * 10) Add control of netCDF header information to the YAML file and this class.
 * 11) Better naming convention for output files?  Currently just keeps a counter, but in EAM the date/time range is used.
 * 12) Better handling of output frequency.  Need to a) gather time info from driver and b) parse that info for frequency units.
 * 13) Expand compatability of IO class to handle Dynamics grids.  This will require a different approach to "set_dofs"
 */

namespace scream
{

class AtmosphereOutput 
{
public:
  using device_type    = DefaultDevice; // may need to template class on this
  using value_type     = Real;          // Right not hard-coded to only Real type variables.  TODO: make this more general?
  using dofs_list_type = KokkosTypes<DefaultDevice>::view_1d<long>;
  using view_type      = typename KokkosTypes<device_type>::template view<value_type*>;

  using Int = scream::Int;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params)
  {
    m_comm   = comm;
    m_params = params;
  }

  // Main Functions
  void init(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm);
  void run(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm, const Real time);
  void finalize();

  // Helper Functions
  void check_status();
  std::map<std::string,Int> get_status() const { return m_status; }

protected:
  // Internal functions
  void add_identifier(const FieldRepository<Real, device_type>& field_repo, const std::string name);
  void register_variables();
  void set_degrees_of_freedom();
  void register_views(const FieldRepository<Real, device_type> &field_repo);
  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;

  // Aaron: You were thinking about making each atmoutput have it's own local field repo where you could store the
  // data needed to handle "averaging" types that are not "instant".  But you keep running into a compilation error.
  std::map<std::string,FieldIdentifier>  m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  dofs_list_type                         m_gids;
  std::map<std::string,view_type>        m_view;

  bool m_is_init = false;

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
inline void AtmosphereOutput::init(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm) 
{
  using namespace scream;
  using namespace scream::scorpio;

  m_status["Init"] += 1;
  m_filename = m_params.get<std::string>("FILENAME")+"_"+std::to_string(m_status["Init"])+".nc";  //TODO: The filename should be treated as a prefix to enable multiple files for the same control.  Like in the case of monthly output with 1 snap/file.
  EKAT_REQUIRE_MSG(!m_is_init,"Error! File " + m_filename + " has already been initialized.\n");

  // Parse the yaml file that controls this output instance.
  m_avg_type        = m_params.get<std::string>("AVERAGING TYPE");
  m_grid_name       = m_params.get<std::string>("GRID");
  auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");

  // Gather data from grid manager:
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with output grid! scorpio.hpp class only supports output on a Physics grid for now.\n");
  m_gids = gm.get_grid(m_grid_name)->get_dofs_gids();
  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  int total_dofs;
  int local_dofs = m_gids.size();
  MPI_Allreduce(&local_dofs, &total_dofs, 1, MPI_INT, MPI_SUM, m_comm.mpi_comm());
  EKAT_REQUIRE_MSG(m_comm.size()<=total_dofs,"Error, PIO interface only allows for the IO comm group size to be less than or equal to the total # of columns in grid.  Consider decreasing size of IO comm group.\n");

  // Create map of fields in this output with the field_identifier in the field repository.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    /* Determine the variable name */
    std::string var_name = var_params.get<std::string>(ekat::strint("field",var_i+1));
    /* Find the <FieldIdentifier,Field> pair that corresponds with this variable name on the appropriate grid */
    add_identifier(field_repo,var_name);
  }

  // Register new netCDF file for output.
  register_outfile(m_filename);

  // Register dimensions with netCDF file.
  for (auto it : m_dims)
  {
    if(it.first == "COL") { it.second=total_dofs;}
    register_dimension(m_filename,it.first,it.first,it.second);
  }
  register_dimension(m_filename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables();
  register_views(field_repo);
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  eam_pio_enddef  (m_filename); 
  m_is_init = true;

} // init
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::run(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm, const Real time) 
{
  using namespace scream;
  using namespace scream::scorpio;

  m_status["Run"] += 1;
  m_status["Avg Count"] += 1;

  // Check to see if output will be written this step
  bool is_write = (m_status["Avg Count"] == m_out_frequency);
  bool avg_init = (m_status["Avg Count"] == 1); // If Avg Count was 0 last step then we need to reset the view.
  // Preamble to writing output this step
  if (is_write) 
  {
    if (!m_is_init) {init(field_repo, gm);} // If m_is_init is false than the previous step was at max_n_steps.  Need a new file.

    pio_update_time(m_filename,time); // Universal scorpio command to set the timelevel for this snap.
    m_status["Snaps"] += 1;           // Update the snap tally, used to determine if a new file is needed.
  }

  // Take care of updating and possibly writing fields.
  for (auto const& f_map : m_fields)
  {
    // Get all the info for this field.
    auto name   = f_map.first;
    auto g_view = field_repo.get_field(f_map.second).get_view();
    auto l_view = m_view.at(name);
    Int  f_len  = f_map.second.get_layout().size();
    // The next two operations do not need to happen if the frequency of output is instantaneous.
    if (m_avg_type != "Instant")
    {
      // If this is the start of a new averaging flag then init the local view.
      if (avg_init) { for (int ii=0; ii<f_len; ++ii) {l_view(ii) = g_view(ii); }
      }
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
    if (is_write) {grid_write_data_array(m_filename,name,m_dofs.at(name),l_view.data());}
  }

  // Finish up any updates to output file and snap counter.
  if (is_write)
  {
    sync_outfile(m_filename);
    // If snaps equals max per file, close this file and set flag to open a new one next write step.
    if (m_status["Snaps"] == m_out_max_steps)
    {
      m_status["Snaps"] = 0;
      eam_pio_closefile(m_filename);
      m_is_init = false;
    }
    // Zero out the Avg Count count now that snap has bee written.
    m_status["Avg Count"] = 0;
  }

} // run
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::finalize() 
{
  using namespace scream;
  using namespace scream::scorpio;
  if (m_is_init) {eam_pio_closefile(m_filename);} // if m_is_init is true then the file was not closed during run step, close it now.
  m_status["Finalize"] += 1;
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::check_status()
{
  printf("IO Status for Rank %5d, File - %.40s: (Init: %2d), (Run: %5d), (Finalize: %2d)\n",m_comm.rank(),m_filename.c_str(),m_status["Init"],m_status["Run"],m_status["Finalize"]);
} // check_status
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::add_identifier(const FieldRepository<Real, device_type> &field_repo, const std::string name)
{
  for (auto myalias=field_repo.begin();myalias!=field_repo.end();++myalias)
  {
    if (myalias->first == name) {
      auto& map_it = myalias->second;
      // map_it is now a map<FieldIdentifier, Field>. Loop over it, and pick all fields defined on the output grid
      // map_it<Field::FieldHeader::FieldIdentifier,Field>
      for (const auto& it : map_it) {
        auto& field_id = it.first;
        if (it.first.get_grid_name()!=m_grid_name) {
          continue;
        }
        // ok, this field is on the correct mesh. add it to pio output
        m_fields.emplace(name,field_id);
        // check to see if all the dims for this field are already set to be registered.
        std::map<std::string,Int>::iterator tag_loc;
        for (int ii=0; ii<field_id.get_layout().rank(); ++ii)
        {
          // check tag against m_dims map.  If not in there, then add it.
          auto& tag = field_id.get_layout().tags()[ii];
          tag_loc = m_dims.find(tag2string(tag));
          if (tag_loc == m_dims.end()) 
          { 
            m_dims.emplace(std::make_pair(tag2string(tag),field_id.get_layout().dim(ii)));
          } else {  
            EKAT_REQUIRE_MSG(m_dims.at(tag2string(tag))==field_id.get_layout().dim(ii),
              "Error! Dimension " + tag2string(tag) + " on field " + name + " has conflicting lengths");
          }
        }
        return;
      }
      // Got this far means we found the field but not on the correct grid.
      EKAT_REQUIRE_MSG(true,"Error! Field " + name + " found in repo, but not on grid " + m_grid_name + ".\n");
    }
  }
  // Got this far means that the field was never found in the repo.
  EKAT_REQUIRE_MSG(true,"Error! Field " + name + " not found in repo.\n");
} // add_identifier
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::register_views(const FieldRepository<Real, device_type> &field_repo)
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    // If the "averaging type" is instant then just need a ptr to the view.
    if (m_avg_type == "Instant")
    {
      m_view.emplace(name,field_repo.get_field(it.second).get_view());
    } else {
    // If anything else then we need a local copy that can be updated.
//      auto view_copy = Kokkos::create_mirror_view( field_repo.get_field(it.second).get_view() );
//      Kokkos::deep_copy(view_copy, field_repo.get_field(it.second).get_view());
      view_type view_copy("",field_repo.get_field(it.second).get_view().extent(0));
      m_view.emplace(name,view_copy);
    }
  }
}
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::register_variables()
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    for (int ii=0;ii<fid.get_layout().rank();++ii) {
      io_decomp_tag += "-" + tag2string(fid.get_layout().tag(ii)); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag2string(fid.get_layout().tag(ii))); // Add dimensions string to vector of dims.
    }
    io_decomp_tag += "-time";  // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger?  Should we register dimension variables (such as ncol and lat/lon) elsewhere in the dimension registration?  These won't have time.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.
    register_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
  }
  // Finish by registering time as a variable.  TODO: Should this really be something registered during the reg. dimensions step? 
  register_variable(m_filename,"time","time",1,{"time"},  PIO_REAL,"t");
} // register_variables
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::set_degrees_of_freedom()
{
  using namespace scream;
  using namespace scream::scorpio;
  // Cycle through all fields and set dof.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
    Int dof_len, n_dim_len;
    // Total number of values represented by this rank for this field is given by the field_layout size.
    dof_len = fid.get_layout().size();
    // For a SCREAM Physics grid, only the total number of columns is decomposed over MPI ranks.  The global id (gid)
    // is stored here as m_gids.  Thus, for this field, the total number of dof's in the other dimensions (i.e. levels)
    // can be found by taking the quotient of dof_len and the length of m_gids.
    n_dim_len = dof_len/m_gids.size();
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global output indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    std::vector<Int> var_dof(dof_len);
    Int dof_it = 0;
    int num_gids = m_gids.size();
    for (int ii=0;ii<num_gids;++ii)
    {
      for (int jj=0;jj<n_dim_len;++jj)
      {
        var_dof[dof_it] =  m_gids(ii)*n_dim_len + jj;
        ++dof_it;
      }
    }
    set_dof(m_filename,name,dof_len,var_dof.data());
    m_dofs.emplace(std::make_pair(name,dof_len));
  }
  // Set degree of freedom for "time"
  set_dof(m_filename,"time",0,0);
  /* TODO: 
   * Adjust DOF to accomodate packing for fields 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
} //namespace scream
#endif // SCREAM_SCORPIO_HPP
