#include "share/io/scorpio_input.hpp"

namespace scream
{

// ====================== IMPLEMENTATION ===================== //
/*
 *  pull_input calls the three typical init, run and finalize routines in order
 *  and makes the input immediately availalbe. Note, init, run and finalize are
 *  also outward facing but do not need those respective sections of the AD, 
 *  i.e. during ad.init, ad.run and ad.finalize.
 */
/* ---------------------------------------------------------- */
void AtmosphereInput::pull_input(const std::string& filename, const std::string& var_name, const std::vector<std::string>& var_dims,
                                        const bool has_columns, const std::vector<int>& dim_lens, const int padding, Real* data)
{
  using namespace scream;
  using namespace scream::scorpio;

  
  auto gids_dev = m_gm->get_grid(m_grid_name)->get_dofs_gids();
  m_gids_host = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids_host,gids_dev);
  // Open the file in PIO
  register_infile(filename);
  // Register the variable information in PIO
  std::string io_decomp_tag = get_io_decomp(var_dims);
  get_variable(filename, var_name, var_name, var_dims.size(), var_dims, PIO_REAL, io_decomp_tag);
  // Determine the degree's of freedom for this variable on this rank, and register with PIO
  int data_length = 1;
  for (auto& dim : dim_lens) { data_length *= dim; }
  std::vector<Int> var_dof = get_var_dof(data_length, has_columns);
  set_dof(filename,var_name,var_dof.size(),var_dof.data());
  set_decomp(filename);
  grid_read_data_array(filename,var_name,dim_lens,data_length,padding,data);
  eam_pio_closefile(filename);  
}
/* ---------------------------------------------------------- */
AtmosphereInput::view_type_host AtmosphereInput::pull_input(const std::string& name)
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  
 *  Overloaded case to deal with just one output and to not put output into field repo.
 *  This is an inefficient way to handle this special case because the same file would need
 *  to be opened and closed multiple times if there is more than one variable.  TODO: Make this
 *  a lot more efficient, most likely will need to overhaul the "pull_input" paradigm.
 */
  using namespace scream;
  using namespace scream::scorpio;
  if (name=="avg_count") {
    view_type_host l_view("",1);
    grid_read_data_array(m_filename,name,m_dofs.at(name),l_view.data());
    return l_view;
  } else {
    auto field = m_field_repo->get_field(name, m_grid_name);
    view_type_host l_view("",field.get_view().extent(0));
    grid_read_data_array(m_filename,name,m_dofs.at(name),l_view.data());
    return l_view;
  }
} 
/* ---------------------------------------------------------- */
void AtmosphereInput::pull_input()
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  */
  using namespace scream;
  using namespace scream::scorpio;
  init();
  for (auto const& name : m_fields)
  {
    auto field = m_field_repo->get_field(name, m_grid_name);
    // Get all the info for this field.
    auto l_dims = field.get_header().get_identifier().get_layout().dims();
    Int padding = field.get_header().get_alloc_properties().get_padding();
    EKAT_REQUIRE_MSG (field.get_header().get_alloc_properties().contiguous(), "Error! Attempting to read in data that is strided for variable: " + name + "\n");
    grid_read_data_array(m_filename,name,l_dims,m_dofs.at(name),padding,field.get_view<Host>().data());
    field.sync_to_dev();
  }
  finalize();
} 
/* ---------------------------------------------------------- */
void AtmosphereInput::init() 
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

  // Check setup against information from grid manager:
  EKAT_REQUIRE_MSG(m_grid_name=="Physics","Error with input grid! scorpio_input.hpp class only supports input on a Physics grid for now.\n");
  auto gids_dev = m_gm->get_grid(m_grid_name)->get_dofs_gids();
  m_gids_host = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids_host,gids_dev); 
  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  int total_dofs;
  int local_dofs = m_gids_host.size();
  MPI_Allreduce(&local_dofs, &total_dofs, 1, MPI_INT, MPI_SUM, m_comm.mpi_comm());
  EKAT_REQUIRE_MSG(m_comm.size()<=total_dofs,"Error, PIO interface only allows for the IO comm group size to be less than or equal to the total # of columns in grid.  Consider decreasing size of IO comm group.\n");

  // Create map of fields in this input with the field_identifier in the field repository.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    /* Determine the variable name */
    std::string var_name = var_params.get<std::string>(ekat::strint("field",var_i+1));
    m_fields.push_back(var_name);
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
void AtmosphereInput::finalize() 
{
/* Cleanup by closing the input file */
  using namespace scream;
  using namespace scream::scorpio;
  eam_pio_closefile(m_filename);
} // finalize
/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
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
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = field.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::vector<std::string> vec_of_dims = get_vec_of_dims(fid.get_layout().tags());
    std::string io_decomp_tag           = get_io_decomp(vec_of_dims);
    get_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire fromt he netCDF the dimensions for each variable.
  }
  if (m_is_rhist) { get_variable(m_filename,"avg_count","avg_count",1,{"cnt"}, PIO_INT, "cnt"); }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<std::string> AtmosphereInput::get_vec_of_dims(const std::vector<FieldTag>& tags)
{
/* Given a set of dimensions in field tags, extract a vector of strings for those dimensions to be used with IO */
  std::vector<std::string> vec_of_dims;
  for (auto& tag_ii : tags)
  {
    vec_of_dims.push_back(e2str(tag_ii));
  }
  std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
  return vec_of_dims;
}
/* ---------------------------------------------------------- */
std::string AtmosphereInput::get_io_decomp(const std::vector<std::string>& vec_of_dims)
{
/* Given a vector of field dimensions, create a unique decomp string to register with I/O/
 * Note: We are hard-coding for only REAL input here.  TODO: would be to allow for other dtypes
 */
  std::string io_decomp_tag = "Real";
  for (auto& tag_ii : vec_of_dims)
  {
    io_decomp_tag += "-" + tag_ii;
  }
  return io_decomp_tag;
}
/* ---------------------------------------------------------- */
void AtmosphereInput::set_degrees_of_freedom()
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
    auto field = m_field_repo->get_field(name, m_grid_name);
    auto& fid  = field.get_header().get_identifier();
    // Given dof_len and n_dim_len it should be possible to create an integer array of "global input indices" for this
    // field and this rank. For every column (i.e. gid) the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    std::vector<Int> var_dof = get_var_dof(fid.get_layout().size(), fid.get_layout().has_tag(FieldTag::Column));
    set_dof(m_filename,name,var_dof.size(),var_dof.data());
    m_dofs.emplace(std::make_pair(name,var_dof.size()));
  }
  if (m_is_rhist)
  {
    Int var_dof[1] = {0};
    set_dof(m_filename,"avg_count",1,var_dof);
    m_dofs.emplace(std::make_pair("avg_count",1));
  }
  /* TODO: 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
std::vector<Int> AtmosphereInput::get_var_dof(const int dof_len, const bool has_cols)
{
  std::vector<Int> var_dof;
  int num_cols;
  // Gather the column degrees of freedom for this variable
  // Total number of values represented by this rank for this field is given by the dof_len.
  // For a SCREAM Physics grid, only the total number of columns is decomposed over MPI ranks.  The global id (gid)
  // is stored here as gids_host. Thus, for this field, the total number of dof's in the other dimensions (i.e. levels)
  // can be found by taking the quotient of dof_len and the length of gids_host.
  if (has_cols)
  {
    num_cols = m_gids_host.size();
  } else {
    // This field is not defined over columns
    // TODO, when we allow for dynamics mesh this check will need to be adjusted for the element tag as well.
    num_cols = 1;
  } 
  Int n_dim_len = dof_len/num_cols;
  // Determine the individual index locations for each degree of freedom.
  for (int ii=0;ii<num_cols;++ii)
  {
    for (int jj=0;jj<n_dim_len;++jj)
    {
      var_dof.push_back(m_gids_host(ii)*n_dim_len + jj);
    }
  }
  return var_dof; 
}
/* ---------------------------------------------------------- */
} // namespace scream
