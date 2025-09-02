#include "share/io/scorpio_input.hpp"

#include "share/field/field_utils.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>

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
AtmosphereInput (const std::string& filename,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::vector<Field>& fields,
                 const bool skip_grid_checks)
{
  // Create param list and field manager on the fly
  ekat::ParameterList params;
  params.set("filename",filename);
  params.set("skip_grid_checks",skip_grid_checks);
  auto& names = params.get<std::vector<std::string>>("field_names",{});

  auto fm = std::make_shared<fm_type>(grid);
  for (auto& f : fields) {
    fm->add_field(f);
    names.push_back(f.name());
  }
  init(params,fm);
}

AtmosphereInput::
AtmosphereInput (const std::string& filename,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::map<std::string,Field>& fields,
                 const bool skip_grid_checks)
{
  // Create param list and field manager on the fly
  ekat::ParameterList params;
  params.set("filename",filename);
  params.set("skip_grid_checks",skip_grid_checks);
  auto& names = params.get<std::vector<std::string>>("field_names",{});

  auto fm = std::make_shared<fm_type>(grid);
  for (const auto& [name, f] : fields) {
    fm->add_field(f);
    names.push_back(name);
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
  finalize();
}

void AtmosphereInput::
set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
  EKAT_REQUIRE_MSG (atm_logger, "Error! Invalid logger pointer.\n");
  m_atm_logger = atm_logger;
}

void AtmosphereInput::
init (const ekat::ParameterList& params,
      const std::shared_ptr<const fm_type>& field_mgr)
{
  EKAT_REQUIRE_MSG (field_mgr->get_grids_manager()->size()==1,
      "Error! AtmosphereInput expects FieldManager defined only on a single grid.\n");
  EKAT_REQUIRE_MSG (not m_fields_inited,
      "Error! Input class was already inited.\n");

  m_params = params;
  m_fields_names = m_params.get<decltype(m_fields_names)>("field_names");
  m_filename = m_params.get<std::string>("filename");

  // Sets the internal field mgr, and possibly sets up the remapper
  set_field_manager(field_mgr);

  // Init scorpio internal structures
  init_scorpio_structures ();
}

void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  // If resetting a field manager we want to check that the layouts of all fields are the same.
  if (m_fm_from_user) {
    for (const auto& [name, f] : m_fm_from_user->get_repo()) {
      auto field_new  = field_mgr->get_field(name);
      // Check Layouts
      auto lay_curr   = f->get_header().get_identifier().get_layout();
      auto lay_new    = field_new.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG(lay_curr==lay_new,"ERROR!! AtmosphereInput::set_field_manager - setting new field manager which has different layout for field " << name <<"\n"
		      << "    Old Layout: " << lay_curr.to_string() << "\n"
		      << "    New Layout: " << lay_new.to_string() << "\n");
    }
  }

  m_fm_from_user = field_mgr;

  // Store grid and fm
  set_grid(field_mgr->get_grid());

  // Init fm_for_scorpio
  m_fm_for_scorpio = std::make_shared<FieldManager>(m_io_grid,RepoState::Closed);
  for (auto const& name : m_fields_names) {
    auto f = m_fm_from_user->get_field(name);
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();

    // If we can alias the field's host view, do it.
    // Otherwise, create a clone.
    bool can_alias = fh.get_parent()==nullptr && fap.get_padding()==0;
    if (can_alias) {
      m_fm_for_scorpio->add_field(f);
    } else {
      // We have padding, or the field is a subfield (or both).
      // Either way, we need a temporary.
      // NOTE: Field::clone honors the max pack size for the allocation,
      //       which would defy one of the reason for cloning it...
      Field copy (f.get_header().get_identifier());
      copy.allocate_view();
      m_fm_for_scorpio->add_field(copy);
    }
  }

  m_fields_inited = true;
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
}

void AtmosphereInput::
reset_filename (const std::string& filename)
{
  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }
  m_params.set("filename",filename);
  m_filename = filename;
  init_scorpio_structures();
}

/* ---------------------------------------------------------- */
void AtmosphereInput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  const bool skip_grid_chk = m_params.get<bool>("skip_grid_checks",false);
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
  m_atm_logger->info("[EAMxx::scorpio_input] Reading variables from file");
  m_atm_logger->info("  file name: " + m_filename);
  m_atm_logger->info("  var names: " + ekat::join(m_fields_names,", "));
  if (time_index!=-1) {
    m_atm_logger->info("  time idx : " + std::to_string(time_index));
  }

  EKAT_REQUIRE_MSG (m_fields_inited and m_scorpio_inited,
      "Error! Internal structures not fully inited yet. Did you forget to call 'init(..)'?\n");

  for (auto const& name : m_fields_names) {

    auto f_scorpio = m_fm_for_scorpio->get_field(name);
    auto f_user    = m_fm_from_user->get_field(name);

    // Read the data
    switch (f_scorpio.data_type()) {
      case DataType::DoubleType:
        scorpio::read_var(m_filename,name,f_scorpio.get_internal_view_data<double,Host>(),time_index);
        break;
      case DataType::FloatType:
        scorpio::read_var(m_filename,name,f_scorpio.get_internal_view_data<float,Host>(),time_index);
        break;
      case DataType::IntType:
        scorpio::read_var(m_filename,name,f_scorpio.get_internal_view_data<int,Host>(),time_index);
        break;
      default:
        EKAT_ERROR_MSG (
            "Error! Unsupported/unrecognized data type while reading field from file.\n"
            " - file name : " + m_filename + "\n"
            " - field name: " + name + "\n");
    }

    f_scorpio.sync_to_dev();
    if (not f_scorpio.is_aliasing(f_user)) {
      f_user.deep_copy(f_scorpio);
    }
  }
  if (m_atm_logger) {
    auto func_finish = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start)/1000.0;
    m_atm_logger->debug("  Done! Elapsed time: " + std::to_string(duration.count()) +" seconds");
  }
}

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize()
{
  if (m_scorpio_inited) {
    scorpio::release_file(m_filename);
  }

  m_fm_from_user   = nullptr;
  m_fm_for_scorpio = nullptr;
  m_io_grid        = nullptr;

  m_fields_inited  = false;
  m_scorpio_inited = false;
}

void AtmosphereInput::init_scorpio_structures()
{
  EKAT_REQUIRE_MSG (m_fields_inited,
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
  for (const auto & [name, f] : m_fm_for_scorpio->get_repo()) {
    const auto& layout = f->get_header().get_identifier().get_layout();

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

  m_scorpio_inited = true;
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
  for (const auto & [name, f] : m_fm_for_scorpio->get_repo()) {
    const auto& layout = f->get_header().get_identifier().get_layout();
    if (layout.has_tag(decomp_tag)) {
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
