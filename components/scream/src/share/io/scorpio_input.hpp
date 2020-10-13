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

namespace scream
{

class AtmosphereInput 
{
public:
  using device_type    = DefaultDevice; // may need to template class on this
  using value_type     = Real;          // Right not hard-coded to only Real type variables.  TODO: make this more general?
  using dofs_list_type = KokkosTypes<DefaultDevice>::view_1d<long>;
  using view_type      = typename KokkosTypes<device_type>::template view<value_type*>;

  using Int = scream::Int;
  using Real = scream::Real;

  virtual ~AtmosphereInput () = default;

  // Constructor
  AtmosphereInput(const ekat::Comm& comm, const ekat::ParameterList& params)
  {
    m_comm   = comm;
    m_params = params;
  }

  // Main Functions
  void pull_input(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm);
  void init(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm);
  void run(const FieldRepository<Real, device_type>& field_repo);
  void finalize();

  // Helper Functions
  void check_status();

protected:
  // Internal functions
  void add_identifier(const FieldRepository<Real, device_type>& field_repo, const std::string name);
  void register_variables();
  void set_degrees_of_freedom();
  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  std::map<std::string,FieldIdentifier>  m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  dofs_list_type                         m_gids;

  bool m_is_init = false;

private:

}; // Class AtmosphereIntput

// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereInput::pull_input(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm)
{
  init(field_repo,gm);
  run(field_repo);
  finalize();
} 
inline void AtmosphereInput::init(const FieldRepository<Real, device_type>& field_repo, const GridsManager& gm) 
{
  using namespace scream;
  using namespace scream::scorpio;

  m_filename = m_params.get<std::string>("FILENAME");

  // Parse the parameters that controls this input instance.
  m_grid_name       = m_params.get<std::string>("GRID");

  // Gather data from grid manager:
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with input grid! scorpio_input.hpp class only supports input on a Physics grid for now.\n");
  m_gids = gm.get_grid(m_grid_name)->get_dofs_gids();
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
    add_identifier(field_repo,var_name);
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
inline void AtmosphereInput::run(const FieldRepository<Real, device_type>& field_repo)
{
  using namespace scream;
  using namespace scream::scorpio;
  using Int = scream::Int;
  using Real = scream::Real;

  // Begin reading all fields and replacing values in field manager
  for (auto const& f_map : m_fields)
  {
    // Get all the info for this field.
    auto name   = f_map.first;
    auto fid    = f_map.second;
    auto l_view_dev = field_repo.get_field(fid).get_view();
    auto l_view_host = Kokkos::create_mirror_view( l_view_dev );
    grid_read_data_array(m_filename,name,m_dofs.at(name),l_view_dev.data());
    Kokkos::deep_copy(l_view_dev,l_view_host);
  }

} // run
/* ---------------------------------------------------------- */
inline void AtmosphereInput::finalize() 
{
  using namespace scream;
  using namespace scream::scorpio;
  eam_pio_closefile(m_filename);
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereInput::add_identifier(const FieldRepository<Real, device_type>& field_repo, const std::string name)
{
  for (auto myalias=field_repo.begin();myalias!=field_repo.end();++myalias)
  {
    if (myalias->first == name) {
      auto& map_it = myalias->second;
      // map_it is now a map<FieldIdentifier, Field>. Loop over it, and pick all fields defined on the input grid
      // map_it<Field::FieldHeader::FieldIdentifier,Field>
      for (const auto& it : map_it) {
        auto& field_id = it.first;
        if (it.first.get_grid_name()!=m_grid_name) {
          continue;
        }
        // ok, this field is on the correct mesh. add it to pio input
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
inline void AtmosphereInput::register_variables()
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
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    get_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire fromt he netCDF the dimensions for each variable.
  }
} // register_variables
/* ---------------------------------------------------------- */
inline void AtmosphereInput::set_degrees_of_freedom()
{
  using namespace scream;
  using namespace scream::scorpio;
  // Cycle through all fields and set dof.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
    bool has_cols = true;
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
      has_cols = false;
    }
    n_dim_len = dof_len/num_cols;
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global input indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    Int var_dof[dof_len];
    Int dof_it = 0;
    for (int ii=0;ii<m_gids.size();++ii)
    {
      for (int jj=0;jj<n_dim_len;++jj)
      {
        var_dof[dof_it] =  m_gids(ii)*n_dim_len + jj;
        ++dof_it;
      }
    }
    set_dof(m_filename,name,dof_len,var_dof);
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
