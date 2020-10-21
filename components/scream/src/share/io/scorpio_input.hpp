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
 *  3. a pointer to the field repository
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

class AtmosphereInput 
{
public:
  using device_type    = DefaultDevice; // may need to template class on this
  using value_type     = Real;          // Right not hard-coded to only Real type variables.  TODO: make this more general?
  using dofs_list_type = KokkosTypes<DefaultDevice>::view_1d<long>;
  using view_type      = typename KokkosTypes<device_type>::template view<value_type*>;

  virtual ~AtmosphereInput () = default;

  // Constructor
  AtmosphereInput(const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const FieldRepository<Real,device_type>>& repo,
                   const std::shared_ptr<const GridsManager>& gm)
  {
    m_comm       = comm;
    m_params     = params;
    m_field_repo = repo;
    m_gm         = gm;
  }

  // Main Functions
  void pull_input();
  void init();
  void run();
  void finalize();

  // Helper Functions
  void check_status();

protected:
  // Internal functions
  void add_identifier(const std::string name);
  void register_variables();
  void set_degrees_of_freedom();
  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  std::shared_ptr<const FieldRepository<Real,device_type>> m_field_repo;
  std::shared_ptr<const GridsManager>                      m_gm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  std::map<std::string,FieldIdentifier>  m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  dofs_list_type                         m_gids;

  bool m_is_init = false;

private:

}; // Class AtmosphereInput

// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereInput::pull_input()
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  */
  init();
  run();
  finalize();
} 
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

  // Gather data from grid manager:
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with input grid! scorpio_input.hpp class only supports input on a Physics grid for now.\n");
  m_gids = m_gm->get_grid(m_grid_name)->get_dofs_gids();
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
    /* Find the <FieldIdentifier,Field> pair that corresponds with this variable name on the appropriate grid */
    add_identifier(var_name);
  }

  // Register new netCDF file for input.
  register_infile(m_filename);

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  set_decomp  (m_filename); 

} // init
/* ---------------------------------------------------------- */
inline void AtmosphereInput::run()
{
/* Replace data in the field repo with data in from the input file for all fields designated as input */
  using namespace scream;
  using namespace scream::scorpio;

  // Begin reading all fields and replacing values in field manager
  for (auto const& f_map : m_fields)
  {
    // Get all the info for this field.
    auto name   = f_map.first;
    auto fid    = f_map.second;
    auto l_view_dev = m_field_repo->get_field(fid).get_view();
    auto l_view_host = Kokkos::create_mirror_view( l_view_dev );
    grid_read_data_array(m_filename,name,m_dofs.at(name),l_view_dev.data());
    Kokkos::deep_copy(l_view_dev,l_view_host);
  }

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
inline void AtmosphereInput::add_identifier(const std::string name)
{
/* 
 * add_identifier routine adds a new field identifier to the local map of fields that will be used by this class for IO.
 * INPUT:
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  auto f = m_field_repo->get_field(name, m_grid_name);
  auto fid = f.get_header().get_identifier();
  m_fields.emplace(name,fid);
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
inline void AtmosphereInput::register_variables()
{
/* Register each variable in IO stream with the SCORPIO interface.  See scream_scorpio_interface.* for details.
 * This is necessary for the SCORPIO routines to be able to lookup variables in the io file with the appropriate
 * degrees of freedom assigned to each core and the correct io decomposition.
 */
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
    for (auto& tag_ii : fid.get_layout().tags()) {
      io_decomp_tag += "-" + tag2string(tag_ii); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag2string(tag_ii)); // Add dimensions string to vector of dims.
    }
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    get_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire fromt he netCDF the dimensions for each variable.
  }
} // register_variables
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
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
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
    std::vector<Int> var_dof(dof_len);
    Int dof_it = 0;
    for (size_t ii=0;ii<m_gids.size();++ii)
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
  /* TODO: 
   * Adjust DOF to accomodate packing for fields 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
} //namespace scream
#endif // SCREAM_SCORPIO_INPUT_HPP
