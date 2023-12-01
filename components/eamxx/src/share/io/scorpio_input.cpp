#include "share/io/scorpio_input.hpp"

#include "share/io/scream_scorpio_interface.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <memory>
#include <numeric>

namespace scream
{

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                 const std::shared_ptr<const fm_type>& field_mgr)
{
  init(params,field_mgr);
}

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::map<std::string,view_1d_host>& host_views_1d,
                 const std::map<std::string,FieldLayout>&  layouts)
{
  init (params,grid,host_views_1d,layouts);
}

AtmosphereInput::
AtmosphereInput (const std::string& filename,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::vector<Field>& fields,
                 const bool skip_grid_checks)
{
  // Create param list and field manager on the fly
  ekat::ParameterList params;
  params.set("Filename",filename);
  params.set("Skip_Grid_Checks",skip_grid_checks);
  auto& names = params.get<std::vector<std::string>>("Field Names",{});

  auto fm = std::make_shared<fm_type>(grid);
  for (auto& f : fields) {
    fm->add_field(f);
    names.push_back(f.name());
  }
  init(params,fm);
}

AtmosphereInput::
~AtmosphereInput ()
{
  // In practice, this should always be true, but since we have a do-nothing default ctor,
  // it is possible to create an instance without ever using it. Since finalize would
  // attempt to close the pio file, we need to call it only if init happened.
  if (m_inited_with_views || m_inited_with_fields) {
    finalize();
  }
}

void AtmosphereInput::
init (const ekat::ParameterList& params,
      const std::shared_ptr<const fm_type>& field_mgr)
{
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited (with user-provided views).\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited (with fields).\n");

  m_params = params;
  m_fields_names = m_params.get<decltype(m_fields_names)>("Field Names");
  m_filename = m_params.get<std::string>("Filename");

  // Sets the internal field mgr, and possibly sets up the remapper
  set_field_manager(field_mgr);

  // Init scorpio internal structures
  init_scorpio_structures ();

  m_inited_with_fields = true;
}

void AtmosphereInput::
init (const ekat::ParameterList& params,
      const std::shared_ptr<const grid_type>& grid,
      const std::map<std::string,view_1d_host>& host_views_1d,
      const std::map<std::string,FieldLayout>&  layouts)
{
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited (with user-provided views).\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited (with fields).\n");

  m_params = params;
  m_filename = m_params.get<std::string>("Filename");

  // Set the grid associated with the input views
  set_grid(grid);

  EKAT_REQUIRE_MSG (host_views_1d.size()==layouts.size(),
      "Error! Input host views and layouts maps has different sizes.\n"
      "       host_views_1d size: " + std::to_string(host_views_1d.size()) + "\n"
      "       layouts size: " + std::to_string(layouts.size()) + "\n");

  m_layouts = layouts;
  m_host_views_1d = host_views_1d;

  // Loop over one of the two maps, store key in m_fields_names,
  // and check that the two maps have the same keys
  for (const auto& it : m_layouts) {
    m_fields_names.push_back(it.first);
    EKAT_REQUIRE_MSG (m_host_views_1d.count(it.first)==1,
        "Error! Input layouts and views maps do not store the same keys.\n"
	"    layout = " + it.first);
  }

  // Init scorpio internal structures
  init_scorpio_structures ();

  m_inited_with_views = true;
}

/* ---------------------------------------------------------- */

void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  // If resetting a field manager we want to check that the layouts of all fields are the same.
  if (m_field_mgr) {
    for (auto felem = m_field_mgr->begin(); felem != m_field_mgr->end(); felem++) {
      auto name = felem->first;
      auto field_curr = m_field_mgr->get_field(name);
      auto field_new  = field_mgr->get_field(name);
      // Check Layouts
      auto lay_curr   = field_curr.get_header().get_identifier().get_layout();
      auto lay_new    = field_new.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG(lay_curr==lay_new,"ERROR!! AtmosphereInput::set_field_manager - setting new field manager which has different layout for field " << name <<"\n"
		      << "    Old Layout: " << to_string(lay_curr) << "\n"
		      << "    New Layout: " << to_string(lay_new) << "\n");
    }
  }

  m_field_mgr = field_mgr;

  // Store grid and fm
  set_grid(m_field_mgr->get_grid());

  // Init fields specs
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

  // The grid is good. Store it.
  m_io_grid = grid;
}

/* ---------------------------------------------------------- */
// Note: The (zero-based) time_index argument provides a way to control which
//       time step to read input from in the file.  If a negative number is
//       provided the routine will read input at the last time level set by
//       running eam_update_timesnap.
void AtmosphereInput::read_variables (const int time_index)
{
  auto func_start = std::chrono::steady_clock::now();
  if (m_atm_logger) {
    m_atm_logger->info("[EAMxx::scorpio_input] Reading variables from file");
    m_atm_logger->info("  file name: " + m_filename);
    m_atm_logger->info("  var names: " + ekat::join(m_fields_names,", "));
    if (time_index!=-1) {
      m_atm_logger->info("  time idx : " + std::to_string(time_index));
    }
  }
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
  auto func_finish = std::chrono::steady_clock::now();
  if (m_atm_logger) {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start)/1000.0;
    m_atm_logger->info("  Done! Elapsed time: " + std::to_string(duration.count()) +" seconds");
  }
} 

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize() 
{
  scorpio::eam_pio_closefile(m_filename);

  m_field_mgr = nullptr;
  m_io_grid   = nullptr;

  m_host_views_1d.clear();
  m_layouts.clear();

  m_inited_with_views = false;
  m_inited_with_fields = false;
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::init_scorpio_structures() 
{
  std::string iotype_str = m_params.get<std::string>("iotype", "default");
  int iotype = scorpio::str2iotype(iotype_str);

  scorpio::register_file(m_filename,scorpio::Read,iotype);

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
    const auto& layout = m_layouts.at(name);
    auto vec_of_dims   = get_vec_of_dims(layout);
    auto io_decomp_tag = get_io_decomp(layout);

    for (size_t  i=0; i<vec_of_dims.size(); ++i) {
      auto partitioned = m_io_grid->get_partitioned_dim_tag()==layout.tags()[i];
      auto dimlen = partitioned ? m_io_grid->get_partitioned_dim_global_size() : layout.dims()[i];
      scorpio::register_dimension(m_filename, vec_of_dims[i], vec_of_dims[i], dimlen, partitioned);
    }

    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to full C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());

    // Register the variable
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire from the netCDF the dimensions for each variable.
    scorpio::register_variable(m_filename, name, name,
                               vec_of_dims, fp_precision, io_decomp_tag);
  }
}

/* ---------------------------------------------------------- */
std::vector<std::string>
AtmosphereInput::get_vec_of_dims(const FieldLayout& layout)
{
  // Given a set of dimensions in field tags, extract a vector of strings
  // for those dimensions to be used with IO
  using namespace ShortFieldTagsNames;
  std::vector<std::string> dims_names;
  dims_names.reserve(layout.rank());
  for (int i=0; i<layout.rank(); ++i) {
    const FieldTag t = layout.tag(i);
    dims_names.push_back(m_io_grid->get_dim_name(layout,i));
    if (t==CMP) {
      dims_names.back() += std::to_string(layout.dim(i));
    }
  }

  return dims_names;
}

/* ---------------------------------------------------------- */
std::string AtmosphereInput::
get_io_decomp(const FieldLayout& layout)
{
  std::string decomp_tag = "dt=real,grid-idx=" + std::to_string(m_io_grid->get_unique_grid_id()) + ",layout=";

  std::vector<int> range(layout.rank());
  std::iota(range.begin(),range.end(),0);
  auto tag_and_dim = [&](int i) {
    return m_io_grid->get_dim_name(layout,i) +
           std::to_string(layout.dim(i));
  };

  decomp_tag += ekat::join (range, tag_and_dim,"-");

  return decomp_tag;
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
  using namespace ShortFieldTagsNames;

  // Precompute this *before* the early return, since it involves collectives.
  // If one rank owns zero cols, and returns prematurely, the others will be left waiting.
  AbstractGrid::gid_type min_gid;
  if (layout.has_tag(COL) or layout.has_tag(EL)) {
    min_gid = m_io_grid->get_global_min_dof_gid();
  }

  // It may be that this MPI ranks owns no chunk of the field
  if (layout.size()==0) {
    return {};
  }

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
  if (layout.has_tag(COL)) {
    const int num_cols = m_io_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      scorpio::offset_t offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else if (layout.has_tag(EL)) {
    auto layout2d = m_io_grid->get_2d_scalar_layout();
    const int num_my_elems = layout2d.dim(0);
    const int ngp = layout2d.dim(1);
    const int num_cols = num_my_elems*ngp*ngp;

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

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
