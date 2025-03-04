#include "share/io/scorpio_input.hpp"

#include "share/io/eamxx_scorpio_interface.hpp"

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
AtmosphereInput (const std::vector<std::string>& fields_names,
                 const std::shared_ptr<const grid_type>& grid)
{
  set_grid(grid);
  m_fields_names = fields_names;
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

  m_inited_with_fields = true;

  // Init scorpio internal structures
  init_scorpio_structures ();
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

  m_inited_with_views = true;

  // Init scorpio internal structures
  init_scorpio_structures ();
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
		      << "    Old Layout: " << lay_curr.to_string() << "\n"
		      << "    New Layout: " << lay_new.to_string() << "\n");
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
    bool can_alias_field_view = fh.get_parent()==nullptr && fap.get_padding()==0;
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

void AtmosphereInput::
set_fields (const std::vector<Field>& fields) {
  auto fm = std::make_shared<fm_type>(m_io_grid);
  m_fields_names.clear();
  for (const auto& f : fields) {
    fm->add_field(f);
    m_fields_names.push_back(f.name());
  }
  set_field_manager(fm);
  m_inited_with_fields = true;
}

void AtmosphereInput::
reset_filename (const std::string& filename)
{
  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }
  m_params.set("Filename",filename);
  m_filename = filename;
  init_scorpio_structures();
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

    scorpio::read_var(m_filename,name,v1d.data(),time_index);

    // If we have a field manager, make sure the data is correctly
    // synced to both host and device views of the field.
    if (m_field_mgr) {

      auto f = m_field_mgr->get_field(name);
      const auto& fh  = f.get_header();
      const auto& fl  = fh.get_identifier().get_layout();
      const auto& fap = fh.get_alloc_properties();

      // Check if the stored 1d view is sharing the data ptr with the field
      const bool can_alias_field_view = fh.get_parent()==nullptr && fap.get_padding()==0;

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
  scorpio::release_file(m_filename);

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
  EKAT_REQUIRE_MSG (m_inited_with_views or m_inited_with_fields,
      "Error! Cannot init scorpio structures until fields/views have been set.\n");

  std::string iotype_str = m_params.get<std::string>("iotype", "default");
  auto iotype = scorpio::str2iotype(iotype_str);

  scorpio::register_file(m_filename,scorpio::Read,iotype);

  // Some input files have the "time" dimension as non-unlimited. This messes up our
  // scorpio interface, which stores a pointer to a "time" dim to be used to read/write
  // slices. This ptr is automatically inited to the unlimited dim in the file. If there is
  // no unlim dim, this ptr remains inited.
  if (not scorpio::has_time_dim(m_filename) and scorpio::has_dim(m_filename,"time")) {
    scorpio::mark_dim_as_time(m_filename,"time");
  }

  // Check variables are in the input file
  for (auto const& name : m_fields_names) {
    const auto& layout = m_layouts.at(name);

    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    auto vec_of_dims   = get_vec_of_dims(layout);

    // Check that the variable is in the file.
    EKAT_REQUIRE_MSG (scorpio::has_var(m_filename,name),
        "Error! Input file does not store a required variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + name + "\n");

    const auto& var = scorpio::get_var(m_filename,name);
    EKAT_REQUIRE_MSG (var.dim_names()==vec_of_dims,
        "Error! Layout mismatch for input file variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + name + "\n"
        " - expected dims : " + ekat::join(vec_of_dims,",") + "\n"
        " - dims from file: " + ekat::join(var.dim_names(),",") + "\n");

    // Check that all dims for this var match the ones on file
    for (int i=0; i<layout.rank(); ++i) {
      const int file_len  = scorpio::get_dimlen(m_filename,vec_of_dims[i]);
      const bool partitioned = m_io_grid->get_partitioned_dim_tag()==layout.tag(i);
      const int eamxx_len = partitioned ? m_io_grid->get_partitioned_dim_global_size()
                                        : layout.dim(i);
      EKAT_REQUIRE_MSG (eamxx_len==file_len,
          "Error! Dimension mismatch for input file variable.\n"
        " - filename : " + m_filename + "\n"
        " - varname  : " + name + "\n"
        " - var dims : " + ekat::join(vec_of_dims,",") + "\n"
        " - dim name : " + vec_of_dims[i] + "\n"
        " - expected extent : " + std::to_string(eamxx_len) + "\n"
        " - extent from file: " + std::to_string(file_len) + "\n");
    }

    // Ensure that we can read the var using Real data type
    scorpio::change_var_dtype (m_filename,name,"real");
  }

  // Set decompositions for the variables
  set_decompositions();
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
    const auto t = layout.tag(i);
    std::string n = m_io_grid->has_special_tag_name(t)
                  ? m_io_grid->get_special_tag_name(t)
                  : layout.names()[i];

    // If t==CMP, and the name stored in the layout is the default ("dim"),
    // we append also the extent, to allow different vector dims in the file
    n += n=="dim" ? std::to_string(layout.dim(i)) : "";

    dims_names.push_back(n);
  }

  return dims_names;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::set_decompositions()
{
  using namespace ShortFieldTagsNames;

  // First, check if any of the vars is indeed partitioned
  const auto decomp_tag  = m_io_grid->get_partitioned_dim_tag();

  bool has_decomposed_layouts = false;
  for (const auto& it : m_layouts) {
    if (it.second.has_tag(decomp_tag)) {
      has_decomposed_layouts = true;
      break;
    }
  }
  if (not has_decomposed_layouts) {
    // If none of the input vars are decomposed on this grid,
    // then there's nothing to do here
    return;
  } 

  // Set the decomposition for the partitioned dimension
  const int local_dim = m_io_grid->get_partitioned_dim_local_size();
  std::string decomp_dim = m_io_grid->has_special_tag_name(decomp_tag)
                         ? m_io_grid->get_special_tag_name(decomp_tag)
                         : e2str(decomp_tag);

  auto gids_f = m_io_grid->get_partitioned_dim_gids();
  auto gids_h = gids_f.get_view<const AbstractGrid::gid_type*,Host>();
  auto min_gid = m_io_grid->get_global_min_partitioned_dim_gid();
  std::vector<scorpio::offset_t> offsets(local_dim);
  for (int idof=0; idof<local_dim; ++idof) {
    offsets[idof] = gids_h[idof] - min_gid;
  }
  scorpio::set_dim_decomp(m_filename,decomp_dim,offsets);
}

} // namespace scream
