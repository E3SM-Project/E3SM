#include "share/io/scorpio_input.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include <memory>
#include <numeric>

namespace scream
{

/* ---------------------------------------------------------- */
AtmosphereInput::
AtmosphereInput (const ekat::Comm& comm,
                 const ekat::ParameterList& params)
 : m_comm (comm)
{
  set_parameters (params);
}

AtmosphereInput::
AtmosphereInput (const ekat::Comm& comm,
                 const ekat::ParameterList& params,
                 const std::shared_ptr<const fm_type>& field_mgr,
                 const std::shared_ptr<const gm_type>& grids_mgr)
 : AtmosphereInput(comm,params)
{
  init (field_mgr,grids_mgr);
}

AtmosphereInput::
AtmosphereInput (const ekat::Comm& comm,
                 const ekat::ParameterList& params,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::map<std::string,view_1d_host>& host_views_1d,
                 const std::map<std::string,FieldLayout>&  layouts)
 : AtmosphereInput(comm,params)
{
  init(grid,host_views_1d,layouts);
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_parameters (const ekat::ParameterList& params) {
  m_filename = params.get<std::string>("Filename");
  m_fields_names = params.get<std::vector<std::string>>("Fields");
}

void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");

  auto fm_grid = field_mgr->get_grid();
  auto unique_grid = fm_grid->get_unique_grid();
  if (fm_grid->is_unique()) {
    // The fm is defined on a unique grid. Store it, then set the grid.
    m_field_mgr = field_mgr;
  } else {
    // The fm is not on a unique grid. We need to create a helper fm,
    // with all the import fields defined on the unique grid,
    // set up scorpio with the helper fm, then create a remapper
    // from the unique grid to the input fm grid.
    EKAT_REQUIRE_MSG (unique_grid,
        "Error! The grid stored in the input field manager is not unique,\n"
        "       and is not storing a unique grid.\n"
        "    grid name: " + fm_grid->name());
    EKAT_REQUIRE_MSG (grids_mgr,
        "Error! The input fields manager is defined on a non-unique grid,\n"
        "       but no grids manager was provided.\n"
        "    grid name: " + fm_grid->name());

    // Create the remapper and register fields from the fid's taken from the input fm
    m_remapper = grids_mgr->create_remapper(unique_grid,fm_grid);
    m_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      auto f = field_mgr->get_field(fname);
      const auto& tgt_fid = f.get_header().get_identifier();

      m_remapper->register_field_from_tgt(tgt_fid);
    }
    m_remapper->registration_ends();

    // Now register the fields in the unique_fm. To generate their
    // field identifiers, use the src fids from the remapper
    auto unique_fm = std::make_shared<fm_type>(unique_grid);
    unique_fm->registration_begins();
    for (int i=0; i<m_remapper->get_num_fields(); ++i) {
      const auto& src_fid = m_remapper->get_src_field_id(i);
      auto f = field_mgr->get_field(src_fid.name());

      const int ps = f.get_header().get_alloc_properties().get_largest_pack_size();
      FieldRequest src_freq(src_fid,ps);
      unique_fm->register_field(src_freq);
    }
    unique_fm->registration_ends();
    m_field_mgr = unique_fm;

    // Finally, bind the src/tgt fields in the remapper
    for (const auto& fname : m_fields_names) {
      auto src = unique_fm->get_field(fname);
      auto tgt = field_mgr->get_field(fname);
      m_remapper->bind_field(src,tgt);
    }

    // This should never fail, but just in case
    EKAT_REQUIRE_MSG (m_remapper->get_num_fields()==m_remapper->get_num_bound_fields(),
        "Error! Something went wrong while building the scorpio input remapper.\n");
  }

  set_grid(unique_grid);
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_grid, "Error! Grid pointer was already set.\n");
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  EKAT_REQUIRE_MSG (grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");

  EKAT_REQUIRE_MSG(m_comm.size()<=grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_grid = grid;
}

/* ---------------------------------------------------------- */
// Note: The time_index argument provides a way to control which
//       time snap to read input from in the file.  If a negative
//       number is provided the routine will read input at the
//       last time level set by running eam_update_timesnap.
void AtmosphereInput::read_variables (const int time_index)
{
  EKAT_REQUIRE_MSG (m_is_inited,
      "Error! The init method has not been called yet.\n");

  for (auto const& name : m_fields_names) {

    // Read the data
    scorpio::grid_read_data_array(m_filename,name,time_index,m_host_views_1d.at(name).data());

    // If we have a field manager, make sure the data is correctly
    // synced to both host and device views of the field.
    if (m_field_mgr) {

      auto f = m_field_mgr->get_field(name);
      const auto& fh  = f.get_header();
      const auto& fl  = fh.get_identifier().get_layout();
      const auto& fap = fh.get_alloc_properties();

      using field_type = decltype(f);
      using RT         = typename field_type::RT;

      // Check if the stored 1d view is sharing the data ptr with the field
      const bool can_alias_field_view = fh.get_parent().expired() && fap.get_padding()==0;

      // If the 1d view is a simple reshape of the field's Host view data,
      // then we're already done. Otherwise, we need to manually copy.
      if (not can_alias_field_view) {
        // Get the host view of the field properly reshaped, and deep copy
        // from temp_view (properly reshaped as well).
        auto rank = fl.rank();
        auto view_1d = m_host_views_1d.at(name);
        switch (rank) {
          case 1:
            {
              // No reshape needed, simply copy
              auto dst = f.get_view<RT*,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                dst(i) = view_1d(i);
              }
              break;
            }
          case 2:
            {
              // Reshape temp_view to a 2d view, then copy
              auto dst = f.get_view<RT**,Host>();
              auto src = view_Nd_host<2>(view_1d.data(),fl.dim(0),fl.dim(1));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  dst(i,j) = src(i,j);
              }}
              break;
            }
          case 3:
            {
              // Reshape temp_view to a 3d view, then copy
              auto dst = f.get_view<RT***,Host>();
              auto src = view_Nd_host<3>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    dst(i,j,k) = src(i,j,k);
              }}}
              break;
            }
          default:
            EKAT_ERROR_MSG ("Error! Unexpected field rank (" + std::to_string(rank) + ").\n");
        }
      }

      // Sync to device
      f.sync_to_dev();
    }
  }

  if (m_remapper) {
    m_remapper->remap(true);
  }
} 

int AtmosphereInput::
read_int_scalar (const std::string& name)
{
  EKAT_REQUIRE_MSG (m_is_inited,
      "Error! The init method has not been called yet.\n");

  return scorpio::get_int_attribute_c2f(m_filename.c_str(),name.c_str());
}

void AtmosphereInput::
init(const std::shared_ptr<const fm_type>& field_mgr,
     const std::shared_ptr<const gm_type>& grids_mgr)
{
  EKAT_REQUIRE_MSG (not m_is_inited,
      "Error! This instance of AtmosphereInput was already inited.\n"
      "       We do not allow re-initialization.\n");

  set_field_manager(field_mgr,grids_mgr);

  for (auto const& name : m_fields_names) {
    auto f = m_field_mgr->get_field(name);
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    const auto& fid = fh.get_identifier();
    const auto& fl  = fid.get_layout();

    // Store tha layout
    m_layouts.emplace(name,fl);

    // If we can alias the field's host view, do it.
    // Otherwise, create a temporary.
    bool can_alias_field_view = fh.get_parent().expired() && fap.get_padding()==0;
    if (can_alias_field_view) {
      auto data = f.get_internal_view_data<Host>();
      m_host_views_1d[name] = view_1d_host(data,fl.size());
    } else {
      // We have padding, or the field is a subfield (or both).
      // Either way, we need a temporary view.
      m_host_views_1d[name] = view_1d_host("",fl.size());
    }
  }

  init_scorpio_structures ();
}

void AtmosphereInput::
init(const std::shared_ptr<const grid_type>& grid,
     const std::map<std::string,view_1d_host>& host_views_1d,
     const std::map<std::string,FieldLayout>&  layouts)
{
  EKAT_REQUIRE_MSG (not m_is_inited,
      "Error! This instance of AtmosphereInput was already inited.\n"
      "       We do not allow re-initialization.\n");

  set_grid(grid);

  // Sanity checks
  EKAT_REQUIRE_MSG (host_views_1d.size()==m_fields_names.size(),
      "Error! Input host views map has the wrong size.\n"
      "       Input size: " + std::to_string(host_views_1d.size()) + "\n"
      "       Expected size: " + std::to_string(m_fields_names.size()) + "\n");
  EKAT_REQUIRE_MSG (layouts.size()==m_fields_names.size(),
      "Error! Input layouts map has the wrong size.\n"
      "       Input size: " + std::to_string(layouts.size()) + "\n"
      "       Expected size: " + std::to_string(m_fields_names.size()) + "\n");

  // Loop over names, rather than just set inputs map in the class.
  // This way, if an expected name is missing, the at(..) method will throw.
  for (auto const& name : m_fields_names) {
    m_layouts.emplace(name,layouts.at(name));
    m_host_views_1d[name] = host_views_1d.at(name);
  }

  init_scorpio_structures ();
}

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize() 
{
  EKAT_REQUIRE_MSG (m_is_inited,
      "Error! The init method has not been called yet.\n");

  scorpio::eam_pio_closefile(m_filename);
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::init_scorpio_structures() 
{
  // This method ensures that the nc file is open, and that all
  // the variables (along with their dimensions) are correctly
  // registered, and sets up scorpio for reading.

  // Register netCDF file for input.
  scorpio::register_file(m_filename,scorpio::Read);

  // Register variables with netCDF file.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  scorpio::set_decomp  (m_filename); 

  m_is_inited = true;
} // init

/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
{
  // Register each variable in IO stream with the SCORPIO interface.
  // This allows SCORPIO to lookup vars in the nc file with the correct
  // dof decomposition across different ranks.

  // Cycle through all fields
  for (auto const& name : m_fields_names) {
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    auto vec_of_dims   = get_vec_of_dims(m_layouts.at(name));
    auto io_decomp_tag = get_io_decomp(vec_of_dims);

    // Register the variable
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire from the netCDF the dimensions for each variable.
    scorpio::get_variable(m_filename, name, name, vec_of_dims.size(),
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
  for (auto it = dims_names.crbegin(); it!=dims_names.crend(); ++it) {
    const auto& dim = *it;
    io_decomp_tag += "-" + dim;
  }

  return io_decomp_tag;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::set_degrees_of_freedom()
{
  // For each field, tell PIO the offset of each DOF to be read.
  // Here, offset is meant in the *global* array in the nc file.
  for (auto const& name : m_fields_names) {
    auto var_dof = get_var_dof_offsets(m_layouts.at(name));
    scorpio::set_dof(m_filename,name,var_dof.size(),var_dof.data());
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

} // namespace scream
