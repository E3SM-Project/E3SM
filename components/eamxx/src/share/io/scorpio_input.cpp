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
 : m_comm   (comm)
 , m_params (params)
{
  m_filename = m_params.get<std::string>("Filename");

  // This ensures that the nc file is open, even if init() doesn't
  // get called. This allows users to read global scalar values from
  // an nc file, by easily creating an AtmosphereInput on the fly.
  scorpio::register_file(m_filename,scorpio::Read);

  // TODO: check that comm is compatible with the pio subsystem comm?
}

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                 const std::shared_ptr<const fm_type>& field_mgr,
                 const std::shared_ptr<const gm_type>& grids_mgr)
 : AtmosphereInput(field_mgr->get_grid()->get_comm(),params)
{
  // Sets the internal field mg, possibly sets up the remapper,
  // and init scorpio internal structures
  init (field_mgr, grids_mgr);
}

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::map<std::string,view_1d_host>& host_views_1d,
                 const std::map<std::string,FieldLayout>&  layouts)
 : AtmosphereInput(grid->get_comm(),params)
{
  // Sets the grid, the host views, and init scorpio internal structures
  init (grid,host_views_1d,layouts);
}

void AtmosphereInput::
init (const std::shared_ptr<const fm_type>& field_mgr,
      const std::shared_ptr<const gm_type>& grids_mgr)
{
  // Set list of fields. Use grid name to potentially find correct sublist inside 'Fields' list.
  set_fields_and_grid_names (field_mgr->get_grid()->aliases());

  // Sets the internal field mgr, and possibly sets up the remapper
  set_field_manager(field_mgr,grids_mgr);

  // Init scorpio internal structures
  init_scorpio_structures ();
}

void AtmosphereInput::
init (const std::shared_ptr<const grid_type>& grid,
      const std::map<std::string,view_1d_host>& host_views_1d,
      const std::map<std::string,FieldLayout>&  layouts)
{
  // Set list of fields. Use grid name to potentially find correct sublist inside 'Fields' list.
  set_fields_and_grid_names (grid->aliases());

  // Set the grid associated with the input views
  set_grid(grid);

  // Set the host views
  set_views(host_views_1d,layouts);

  // Init scorpio internal structures
  init_scorpio_structures ();
}

/* ---------------------------------------------------------- */

void AtmosphereInput::
set_fields_and_grid_names (const std::vector<std::string>& grid_aliases) {
  // The user might just want to read some global attributes (no fields),
  // so get the list of fields names only if present.
  using vos_t = std::vector<std::string>;
  if (m_params.isParameter("Field Names")) {
    m_fields_names = m_params.get<vos_t>("Field Names");
    if (m_params.isParameter("IO Grid Name")) {
      m_io_grid_name = m_params.get<std::string>("IO Grid Name");
    }
  } else {
    for (const auto& grid_name : grid_aliases) {
      if (m_params.isSublist("Fields") && grid_name!="") {
        const auto& pl = m_params.sublist("Fields").sublist(grid_name);
        m_fields_names = pl.get<vos_t>("Field Names");
        if (pl.isParameter("IO Grid Name")) {
          m_io_grid_name = pl.get<std::string>("IO Grid Name");
        }
      }
      break;
    }
  }
}

void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited with user-provided views.\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited with fields.\n");

  m_field_mgr = field_mgr;

  std::shared_ptr<const grid_type> fm_grid, io_grid;
  io_grid = fm_grid = m_field_mgr->get_grid();

  if (m_io_grid_name!="" && m_io_grid_name!=fm_grid->name()) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    io_grid = grids_mgr->get_grid(m_io_grid_name);
    m_remapper = grids_mgr->create_remapper(io_grid,fm_grid);

    // Register all input fields in the remapper.
    m_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      auto f = m_field_mgr->get_field(fname);
      const auto& tgt_fid = f.get_header().get_identifier();
      EKAT_REQUIRE_MSG(tgt_fid.data_type()==DataType::RealType,
          "Error! I/O supports only Real data, for now.\n");
      m_remapper->register_field_from_tgt(tgt_fid);
    }
    m_remapper->registration_ends();

    // Now create a new FM on io grid, and create copies of input fields from FM.
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (int i=0; i<m_remapper->get_num_fields(); ++i) {
      const auto& src_fid = m_remapper->get_src_field_id(i);
      io_fm->register_field(FieldRequest(src_fid));
    }
    io_fm->registration_ends();

    // Now that fields have been allocated on the io grid, we can bind them in the remapper
    for (const auto& fname : m_fields_names) {
      auto src = io_fm->get_field(fname);
      auto tgt = m_field_mgr->get_field(fname);
      m_remapper->bind_field(src,tgt);
    }

    // This should never fail, but just in case
    EKAT_REQUIRE_MSG (m_remapper->get_num_fields()==m_remapper->get_num_bound_fields(),
        "Error! Something went wrong while building the scorpio input remapper.\n");

    // Reset field mgr
    m_field_mgr = io_fm;
  }

  // Store grid and fm
  set_grid(io_grid);

  // Init fields specs
  register_fields_specs();

  m_inited_with_fields = true;
}

void AtmosphereInput::
register_fields_specs() {
  for (auto const& name : m_fields_names) {
    auto f = m_field_mgr->get_field(name);
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    const auto& fid = fh.get_identifier();
    const auto& fl  = fid.get_layout();

    // Store the layout
    m_layouts.emplace(name,fl);

    // If we can alias the field's host view, do it.
    // Otherwise, create a temporary.
    bool can_alias_field_view = fh.get_parent().expired() && fap.get_padding()==0;
    if (can_alias_field_view) {
      auto data = f.get_internal_view_data<Real,Host>();
      m_host_views_1d[name] = view_1d_host(data,fl.size());
    } else {
      // We have padding, or the field is a subfield (or both).
      // Either way, we need a temporary view.
      m_host_views_1d[name] = view_1d_host("",fl.size());
    }
  }
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_io_grid, "Error! Grid pointer was already set.\n");
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  const bool skip_grid_chk = m_params.get<bool>("Skip_Grid_Checks",false);
  if (!skip_grid_chk) {
    EKAT_REQUIRE_MSG (grid->is_unique(),
        "Error! I/O only supports grids which are 'unique', meaning that the\n"
        "       map dof_gid->proc_id is well defined.\n");
    EKAT_REQUIRE_MSG (
        (grid->get_global_max_dof_gid()-grid->get_global_min_dof_gid()+1)==grid->get_num_global_dofs(),
        "Error! IO requires DOF gids to (globally)  be in interval [gid_0,gid_0+num_global_dofs).\n"
        "   - global min GID : " + std::to_string(grid->get_global_min_dof_gid()) + "\n"
        "   - global max GID : " + std::to_string(grid->get_global_max_dof_gid()) + "\n"
        "   - num global dofs: " + std::to_string(grid->get_num_global_dofs()) + "\n");
  }

  EKAT_REQUIRE_MSG(m_comm.size()<=grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_io_grid = grid;

  // Reset the comm
  m_comm = m_io_grid->get_comm();
}

/* ---------------------------------------------------------- */
// Note: The (zero-based) time_index argument provides a way to control which
//       time step to read input from in the file.  If a negative number is
//       provided the routine will read input at the last time level set by
//       running eam_update_timesnap.
void AtmosphereInput::read_variables (const int time_index)
{
  EKAT_REQUIRE_MSG (m_inited_with_views || m_inited_with_fields,
      "Error! Scorpio structures not inited yet. Did you forget to call 'init(..)'?\n");

  for (auto const& name : m_fields_names) {

    // Read the data
    auto v1d = m_host_views_1d.at(name);
    scorpio::grid_read_data_array(m_filename,name,time_index,v1d.data(),v1d.size());

    // If we have a field manager, make sure the data is correctly
    // synced to both host and device views of the field.
    if (m_field_mgr) {

      auto f = m_field_mgr->get_field(name);
      const auto& fh  = f.get_header();
      const auto& fl  = fh.get_identifier().get_layout();
      const auto& fap = fh.get_alloc_properties();

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
              auto dst = f.get_view<Real*,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                dst(i) = view_1d(i);
              }
              break;
            }
          case 2:
            {
              // Reshape temp_view to a 2d view, then copy
              auto dst = f.get_view<Real**,Host>();
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
              auto dst = f.get_view<Real***,Host>();
              auto src = view_Nd_host<3>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    dst(i,j,k) = src(i,j,k);
              }}}
              break;
            }
          case 4:
            {
              // Reshape temp_view to a 4d view, then copy
              auto dst = f.get_view<Real****,Host>();
              auto src = view_Nd_host<4>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      dst(i,j,k,l) = src(i,j,k,l);
              }}}}
              break;
            }
          case 5:
            {
              // Reshape temp_view to a 5d view, then copy
              auto dst = f.get_view<Real*****,Host>();
              auto src = view_Nd_host<5>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3),fl.dim(4));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      for (int m=0; m<fl.dim(4); ++m) {
                        dst(i,j,k,l,m) = src(i,j,k,l,m);
              }}}}}
              break;
            }
          case 6:
            {
              // Reshape temp_view to a 6d view, then copy
              auto dst = f.get_view<Real******,Host>();
              auto src = view_Nd_host<6>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3),fl.dim(4),fl.dim(5));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      for (int m=0; m<fl.dim(4); ++m) {
                        for (int n=0; n<fl.dim(5); ++n) {
                          dst(i,j,k,l,m,n) = src(i,j,k,l,m,n);
              }}}}}}
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
  return scorpio::get_int_attribute_c2f(m_filename.c_str(),name.c_str());
}

void AtmosphereInput::
set_views (const std::map<std::string,view_1d_host>& host_views_1d,
           const std::map<std::string,FieldLayout>&  layouts)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited with user-provided views.\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited with fields.\n");
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

  m_inited_with_views = true;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize() 
{
  scorpio::eam_pio_closefile(m_filename);

  m_field_mgr = nullptr;
  m_io_grid   = nullptr;
  m_remapper  = nullptr;

  m_host_views_1d.clear();
  m_layouts.clear();

  m_inited_with_views = false;
  m_inited_with_fields = false;
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::init_scorpio_structures() 
{
  // Register variables with netCDF file.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  scorpio::set_decomp  (m_filename); 
}

/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
{
  // Register each variable in IO stream with the SCORPIO interface.
  // This allows SCORPIO to lookup vars in the nc file with the correct
  // dof decomposition across different ranks.

  // Cycle through all fields
  const auto& fp_precision = "real";
  for (auto const& name : m_fields_names) {
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    auto vec_of_dims   = get_vec_of_dims(m_layouts.at(name));
    auto io_decomp_tag = get_io_decomp(m_layouts.at(name));

    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to full C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());

    // Register the variable
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire from the netCDF the dimensions for each variable.
    scorpio::get_variable(m_filename, name, name,
                          vec_of_dims, fp_precision, io_decomp_tag);
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

  return dims_names;
}

/* ---------------------------------------------------------- */
std::string AtmosphereInput::
get_io_decomp(const FieldLayout& layout)
{
  // Given a vector of dimensions names, create a unique decomp string to register with I/O
  // Note: We are hard-coding for only REAL input here.
  // TODO: would be to allow for other dtypes
  std::string io_decomp_tag = (std::string("Real-") + m_io_grid->name() + "-" +
                               std::to_string(m_io_grid->get_num_global_dofs()));
  auto dims_names = get_vec_of_dims(layout);
  for (size_t i=0; i<dims_names.size(); ++i) {
    io_decomp_tag += "-" + dims_names[i];
    // If tag==CMP, we already attached the length to the tag name
    if (layout.tag(i)!=ShortFieldTagsNames::CMP) {
      io_decomp_tag += "_" + std::to_string(layout.dim(i));
    }
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
std::vector<scorpio::offset_t>
AtmosphereInput::get_var_dof_offsets(const FieldLayout& layout)
{
  std::vector<scorpio::offset_t> var_dof(layout.size());

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // Since we order the global array based on dof gid, and we *assume* (we actually
  // check this during set_grid) that the grid global gids are in the interval
  // [gid_0, gid_0+num_global_dofs), the offset is simply given by
  // (dof_gid-gid_0)*column_size (for partitioned arrays).
  // NOTE: a "dof" in the grid object is not the same as a "dof" in scorpio.
  //       For a SEGrid 3d vector field with (MPI local) layout (nelem,2,np,np,nlev),
  //       scorpio sees nelem*2*np*np*nlev dofs, while the SE grid sees nelem*np*np dofs.
  //       All we need to do in this routine is to compute the offset of all the entries
  //       of the MPI-local array w.r.t. the global array. So long as the offsets are in
  //       the same order as the corresponding entry in the data to be read/written, we're good.
  auto dofs_h = m_io_grid->get_dofs_gids().get_view<const AbstractGrid::gid_type*,Host>();
  if (layout.has_tag(ShortFieldTagsNames::COL)) {
    const int num_cols = m_io_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      scorpio::offset_t offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else if (layout.has_tag(ShortFieldTagsNames::EL)) {
    auto layout2d = m_io_grid->get_2d_scalar_layout();
    const int num_my_elems = layout2d.dim(0);
    const int ngp = layout2d.dim(1);
    const int num_cols = num_my_elems*ngp*ngp;

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int ie=0,icol=0; ie<num_my_elems; ++ie) {
      for (int igp=0; igp<ngp; ++igp) {
        for (int jgp=0; jgp<ngp; ++jgp,++icol) {
          // Get chunk of var_dof to fill
          auto start = var_dof.begin()+icol*col_size;
          auto end   = start+col_size;

          // Compute start of the column offset, then fill column adding 1 to each entry
          auto gid = dofs_h(icol);
          auto offset = (gid-min_gid)*col_size;
          std::iota(start,end,offset);
    }}}
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}

} // namespace scream
