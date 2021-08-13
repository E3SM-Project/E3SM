#include "share/io/scorpio_input.hpp"

#include "ekat/ekat_parameter_list.hpp"

#include <numeric>

namespace scream
{

/* ---------------------------------------------------------- */
AtmosphereInput::
AtmosphereInput (const ekat::Comm& comm,
                 const ekat::ParameterList& params,
                 const std::shared_ptr<const fm_type>& field_mgr)
 : m_comm (comm)
{
  set_parameters (params);
  set_field_manager (field_mgr);

  init ();
}

/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_parameters (const ekat::ParameterList& params) {
  m_params = params;

  m_filename = m_params.get<std::string>("FILENAME");
  m_fields_names = m_params.get<std::vector<std::string>>("FIELDS");

  // If this is a history restart type of read, make sure its noted
  m_is_history_restart = m_params.get("History Restart", false);
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");

  // The fm is good. Store it, then set the grid.
  m_field_mgr = field_mgr;
  set_grid(m_field_mgr->get_grid());
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_grid, "Error! Grid pointer was already set.\n");
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  EKAT_REQUIRE_MSG (grid->name()=="Physics" || grid->name()=="Physics GLL",
      "Error! I/O only supports output on a Physics or Physics GLL grid.\n");

  EKAT_REQUIRE_MSG(m_comm.size()<=grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_grid = grid;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::read_variables ()
{
  using namespace scream::scorpio;

  for (auto const& name : m_fields_names) {
    // Retrieve field, and some of its metadata
    auto field = m_field_mgr->get_field(name);
    const auto& fh  = field.get_header();
    const auto& fl  = fh.get_identifier().get_layout();
    const auto& fap = fh.get_alloc_properties();

    // The field might be a subfield or have padding, in which case its view
    // is not a contiguous array in memory.
    // Strategy: create a temp contiguous view, and use it for reading from file.
    // Then, deep copy back to the input view. If the field does *not* have
    // a parent and is *not* padded, the temp view can store the same memory
    // pointer as the Host view in the Field.
    using field_type = decltype(field);
    using RT         = typename field_type::RT;

    // Use a 1d view of correct (i.e., physical) size for scorpio reading
    ekat::Unmanaged<view_1d_host> temp_view;
    view_1d_host scratch_view;
    bool needs_scratch_view = not (fh.get_parent().expired() && fap.get_padding()==0);
    if (needs_scratch_view) {
      scratch_view = decltype(scratch_view)("",fl.size());
      temp_view = decltype(temp_view)(scratch_view.data(),fl.size());
    } else {
      temp_view = decltype(temp_view)(field.get_internal_view_data<Host>(),fl.size());
    }

    // Read the data
    grid_read_data_array(m_filename,name,temp_view.data());

    // If temp_view is a simple reshape of the field's Host view data, then we're alreaddy done.
    // Otherwise, we need to do a deep copy.
    if (needs_scratch_view) {
      // Get the host view of the field properly reshaped, and deep copy
      // from temp_view (properly reshaped as well).
      auto rank = fl.rank();
      switch (rank) {
        case 1:
          {
            // No reshape needed, simply copy
            auto dst = field.get_view<RT*,Host>();
            for (int i=0; i<fl.dim(0); ++i) {
              dst(i) = temp_view(i);
            }
            break;
          }
        case 2:
          {
            // Reshape temp_view to a 2d view, then copy
            auto dst = field.get_view<RT**,Host>();
            auto src = view_ND_host<2>(temp_view.data(),fl.dim(0),fl.dim(1));
            Kokkos::Impl::ViewRemap<decltype(dst),decltype(src),HostDevice,2>(dst,src);
            break;
          }
        case 3:
          {
            // Reshape temp_view to a 3d view, then copy
            auto dst = field.get_view<RT***,Host>();
            auto src = view_ND_host<3>(temp_view.data(),fl.dim(0),fl.dim(1),fl.dim(2));
            Kokkos::Impl::ViewRemap<decltype(dst),decltype(src),HostDevice,3>(dst,src);
            break;
          }
        default:
          EKAT_ERROR_MSG ("Error! Unexpected field rank (" + std::to_string(rank) + ").\n");
      }
    }

    // Sync to device
    field.sync_to_dev();
  }
} 

int AtmosphereInput::
read_int_scalar (const std::string& name)
{
  return scorpio::get_int_attribute_c2f(m_filename.c_str(),name.c_str());
}

/* ---------------------------------------------------------- */
void AtmosphereInput::init() 
{
  // This method ensures that the nc file is open, and that all
  // the variables (along with their dimensions) are correctly
  // registered, and sets up scorpio for reading.

  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // Sanity check: if a grid was specified in the input file,
  // it must match the one stored.
  if (m_params.isParameter("GRID")) {
    EKAT_REQUIRE_MSG (m_params.get<std::string>("GRID")==m_grid->name(),
        "Error! Input grid name in the parameter list does not match the name of the input grid.\n");
  }

  // Register netCDF file for input.
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
  scorpio::eam_pio_closefile(m_filename);
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
{
  // Register each variable in IO stream with the SCORPIO interface.
  // This allows SCORPIO to lookup vars in the nc file with the correct
  // dof decomposition across different ranks.

  using namespace scorpio;

  // Cycle through all fields
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();

    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    auto vec_of_dims   = get_vec_of_dims(fid.get_layout());
    auto io_decomp_tag = get_io_decomp(vec_of_dims);

    // Register the variable
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire from the netCDF the dimensions for each variable.
    get_variable(m_filename, name, name, vec_of_dims.size(),
                 vec_of_dims, PIO_REAL, io_decomp_tag);
  }
}

/* ---------------------------------------------------------- */
std::vector<std::string>
AtmosphereInput::get_vec_of_dims(const FieldLayout& layout)
{
  // Given a set of dimensions in field tags, extract a vector of strings
  // for those dimensions to be used with IO
  std::vector<std::string> dims_names(layout.rank());
  for (int i=0; i<layout.rank(); ++i) {
    dims_names[i] = scorpio::get_nc_tag_name(layout.tag(i),layout.dim(i));
  }

  // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
  // may need to delete this line when switching to full C++/C implementation.
  std::reverse(dims_names.begin(),dims_names.end());

  return dims_names;
}

/* ---------------------------------------------------------- */
std::string AtmosphereInput::get_io_decomp(const std::vector<std::string>& dims_names)
{
  // Given a vector of dimensions names, create a unique decomp string to register with I/O
  // Note: We are hard-coding for only REAL input here.
  // TODO: would be to allow for other dtypes
  std::string io_decomp_tag = "Real";
  for (const auto& dim : dims_names) {
    io_decomp_tag += "-" + dim;
  }

  return io_decomp_tag;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::set_degrees_of_freedom()
{
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // For each field, tell PIO the offset of each DOF to be read.
  // Here, offset is meant in the *global* array in the nc file.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();

    auto var_dof = get_var_dof_offsets(fid.get_layout());
    set_dof(m_filename,name,var_dof.size(),var_dof.data());
  }
} // set_degrees_of_freedom

/* ---------------------------------------------------------- */
std::vector<int> AtmosphereInput::
get_var_dof_offsets(const FieldLayout& layout)
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
    Int col_size = layout.size() / num_cols;

    // Compute the number of columns owned by all previous ranks.
    int offset = 0;
    m_comm.scan_sum(&num_cols,&offset,1);

    // Compute offsets of all my dofs
    std::iota(var_dof.begin(), var_dof.end(), offset*col_size);
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}

} // namespace scream
