#ifndef SCREAM_SCORPIO_INPUT_HPP
#define SCREAM_SCORPIO_INPUT_HPP

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
/*  The AtmosphereInput class handles all input streams to SCREAM.
 *  It is important to note that there does not exist an InputManager,
 *  like in the case of output.  So all input streams have to be managed
 *  directly by the process that requires it.
 *
 *  The typical input call will be to the outward facing routine 'pull_input'.
 *
 *  pullinput calls the three typical init, run and finalize routines in order
 *  and makes the input immediately availalbe. Note, init, run and finalize are
 *  also outward facing but do not need those respective sections of the AD, 
 *  i.e. during ad.init, ad.run and ad.finalize.
 *
 *  Currently, input can only handle single timesnap input files.  In other words
 *  files that will be opened, read and closed within the same timestep and only
 *  store one timesnap of data.
 *  TODO: Develop infrastructure for reading input from files with multiple timesnaps
 *  and over multiple timesteps of the simulation.
 *
 *  Note: the init, run and finalize are separate routines that are outward facing in
 *  this class to facilitate cases where reading input over some number of simulation
 *  timesteps is possible.
 *
 *  At construction time All input instances require
 *  1. an EKAT comm group and
 *  2. a EKAT parameter list
 *  3. a shared pointer to the field repository
 *  4. a shared pointer to the grids manager
 *  The parameter list contains all of the essential data regarding
 *  the input file name and the set of variables to be read.  The parameter list can be
 *  created localling in the process needing input, see src/share/io/tests/ for examples of
 *  setting up the input parameter list.
 *
 *  A typical input parameter list looks like:
 *  -----
 *  Input Parameters
 *    FILENAME: STRING
 *    GRID: STRING
 *    FIELDS
 *      Number of Fields: INT
 *      field_1: STRING
 *      ...
 *      field_N: STRING
 *  -----
 *  where,
 *  FILENAME: is a string value of the name of the input file to be read.
 *  GRID: is a string of the grid the input file is written on, currently only "Physics" is supported.
 *  FIELDS: designation of a sublist, so empty here
 *    Number of Fields: is an integer value>0 telling the number of fields
 *    field_x: is the xth field variable name.  Should match the name in the file and the name in the field repo.  TODO: add a rename option if variable names differ in file and field repo.
 *
 * Usage:
 * 1. Construct an instance of the AtmosphereInput class:
 *    AtmosphereInput in_object(comm,params);
 * 2. Use pull input to gather the desired data:
 *    in_object.pull_input(*grid_manager)
 *
 * The AtmosphereInput class will replace all fields in the field_repo that are part of the input with
 * data read from the input file.
 *    
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */
namespace scream
{

  using field_type      = Field<Real>;
  using header_type     = typename field_type::header_type;
  using identifier_type = typename header_type::identifier_type;

class AtmosphereInput 
{
public:
  using dofs_list_type = KokkosTypes<HostDevice>::view_1d<long>;
  using view_type      = typename KokkosTypes<HostDevice>::view_1d<Real>;

  virtual ~AtmosphereInput () = default;

  // Constructor
  AtmosphereInput(const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const FieldRepository<Real>>& repo,
                  const std::shared_ptr<const GridsManager>& gm)
  {
    m_comm       = comm;
    m_params     = params;
    m_field_repo = repo;
    m_gm         = gm;
  }

  // Main Functions
  void pull_input();
  void pull_input(const std::string name, view_type view_out);
  void init();
  view_type run(const std::string name);
  void finalize();

  // Helper Functions
  void check_status();

protected:
  // Internal functions
  void register_dimensions(const std::string name);
  void register_variables();
  void set_degrees_of_freedom();
  void register_views();
  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  std::shared_ptr<const FieldRepository<Real>> m_field_repo;
  std::shared_ptr<const GridsManager>                      m_gm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  std::vector<std::string>               m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  typename dofs_list_type::HostMirror    m_gids;
  std::map<std::string,view_type>        m_view_local;

  bool m_is_init  = false;
  bool m_is_rhist = false;

private:

}; // Class AtmosphereInput

// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereInput::pull_input(const std::string name, view_type view_out)
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  
 *  Overloaded case to deal with just one output and to not put output into field repo.
 *  This is an inefficient way to handle this special case because the same file would need
 *  to be opened and closed multiple times if there is more than one variable.  TODO: Make this
 *  a lot more efficient, most likely will need to overhaul the "pull_input" paradigm.
 */
  auto l_view = run(name);
  Kokkos::deep_copy(view_out,l_view);
} 
/* ---------------------------------------------------------- */
inline void AtmosphereInput::pull_input()
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  */
  init();
  for (auto const& name : m_fields)
  {
    auto fmap = m_field_repo->get_field(name, m_grid_name);
    // Get all the info for this field.
    auto l_view_dev = fmap.get_view();
    auto l_view = run(name);
    Kokkos::deep_copy(l_view_dev,l_view);
  }
  finalize();
} 
/* ---------------------------------------------------------- */
inline void AtmosphereInput::init() 
{
/* 
 * Call all of the necessary SCORPIO routines to open the file, gather the variables and dimensions and set the degrees-of-freedom
 * for reading.
 * Organize local structures responsible for writing over field repo data.
 */
  using namespace scream;
  using namespace scream::scorpio;

  m_filename = m_params.get<std::string>("FILENAME");

  // Parse the parameters that controls this input instance.
  m_grid_name       = m_params.get<std::string>("GRID");

  // If this is RHIST type make sure its noted
  if (m_params.isParameter("RHIST"))
  {
    m_is_rhist = m_params.get<bool>("RHIST");
  }

  // Gather data from grid manager:
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with input grid! scorpio_input.hpp class only supports input on a Physics grid for now.\n");
  auto gids_dev = m_gm->get_grid(m_grid_name)->get_dofs_gids();
  m_gids = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids,gids_dev); 
  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  int total_dofs;
  int local_dofs = m_gids.size();
  MPI_Allreduce(&local_dofs, &total_dofs, 1, MPI_INT, MPI_SUM, m_comm.mpi_comm());
  EKAT_REQUIRE_MSG(m_comm.size()<=total_dofs,"Error, PIO interface only allows for the IO comm group size to be less than or equal to the total # of columns in grid.  Consider decreasing size of IO comm group.\n");

  // Create map of fields in this input with the field_identifier in the field repository.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    /* Determine the variable name */
    std::string var_name = var_params.get<std::string>(ekat::strint("field",var_i+1));
    m_fields.push_back(var_name);
    /* Check that all dimensions for this variable are set to be registered */
    register_dimensions(var_name);
  }

  // Now that the fields have been gathered register the local views which will be used to determine output data to be written.
  register_views();

  // Register new netCDF file for input.
  register_infile(m_filename);

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  set_decomp  (m_filename); 

} // init
/* ---------------------------------------------------------- */
inline AtmosphereInput::view_type AtmosphereInput::run(const std::string name)
{
/* Load data from the input file */
  using namespace scream;
  using namespace scream::scorpio;

    view_type l_view_host = m_view_local.at(name);
    grid_read_data_array(m_filename,name,m_dofs.at(name),l_view_host.data());
    return l_view_host;

} // run
/* ---------------------------------------------------------- */
inline void AtmosphereInput::finalize() 
{
/* Cleanup by closing the input file */
  using namespace scream;
  using namespace scream::scorpio;
  eam_pio_closefile(m_filename);
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereInput::register_dimensions(const std::string name)
{
/* 
 * Checks that the dimensions associated with a specific variable will be registered with IO file.
 * INPUT:
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  auto f = m_field_repo->get_field(name, m_grid_name);
  auto fid = f.get_header().get_identifier();
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
} // register_dimensions
/* ---------------------------------------------------------- */
inline void AtmosphereInput::register_variables()
{
/* Register each variable in IO stream with the SCORPIO interface.  See scream_scorpio_interface.* for details.
 * This is necessary for the SCORPIO routines to be able to lookup variables in the io file with the appropriate
 * degrees of freedom assigned to each core and the correct io decomposition.
 */
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields)
  {
    auto fmap = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = fmap.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    for (auto& tag_ii : fid.get_layout().tags()) {
      io_decomp_tag += "-" + tag2string(tag_ii); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag2string(tag_ii)); // Add dimensions string to vector of dims.
    }
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    get_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire fromt he netCDF the dimensions for each variable.
  }
  if (m_is_rhist) { get_variable(m_filename,"avg_count","avg_count",1,{"cnt"}, PIO_INT, "cnt"); }
} // register_variables
/* ---------------------------------------------------------- */
inline void AtmosphereInput::register_views()
{
  using namespace scream;
  using namespace scream::scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields)
  {
    auto fmap = m_field_repo->get_field(name, m_grid_name);
    // If the "averaging type" is instant then just need a ptr to the view.
    auto view_d = fmap.get_view();
    // Create a local copy of view to be stored by input stream.
    auto view_copy_h = Kokkos::create_mirror_view( view_d );
    view_type view_copy("",view_copy_h.extent(0));
    Kokkos::deep_copy(view_copy, view_d);
    m_view_local.emplace(name,view_copy);
  }
  if (m_is_rhist)
  {
    view_type l_view("",1);
    m_view_local.emplace("avg_count",l_view);
  }
}
/* ---------------------------------------------------------- */
inline void AtmosphereInput::set_degrees_of_freedom()
{
/* 
 * Use information from the grids manager to determine which indices for each field are associated
 * with this computational core.
 */
  using namespace scream;
  using namespace scream::scorpio;
  // Cycle through all fields and set dof.
  for (auto const& name : m_fields)
  {
    auto fmap = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = fmap.get_header().get_identifier();
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
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global input indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    std::vector<Int> var_dof;
    for (size_t ii=0;ii<m_gids.size();++ii)
    {
      for (int jj=0;jj<n_dim_len;++jj)
      {
        var_dof.push_back(m_gids(ii)*n_dim_len + jj);
      }
    }
    set_dof(m_filename,name,dof_len,var_dof.data());
    m_dofs.emplace(std::make_pair(name,dof_len));
  }
  if (m_is_rhist)
  {
    Int var_dof[1] = {0};
    set_dof(m_filename,"avg_count",1,var_dof);
    m_dofs.emplace(std::make_pair("avg_count",1));
  }
  /* TODO: 
   * Adjust DOF to accomodate packing for fields 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
} //namespace scream
#endif // SCREAM_SCORPIO_INPUT_HPP
