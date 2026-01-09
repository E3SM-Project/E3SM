#include "share/data_managers/field_reader.hpp"

#include "share/field/field_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>

#include <memory>
#include <numeric>

namespace scream
{

FieldReader::
FieldReader (const std::string& filename,
             const std::shared_ptr<const AbstractGrid>& grid,
             const std::vector<Field>& fields,
             const std::string& iotype)
{
  m_filename = filename;

  set_grid(grid);

  set_fields(fields);

  init_scorpio_structures(iotype);
}

FieldReader::
FieldReader (const std::shared_ptr<const AbstractGrid>& grid)
{
  set_grid(grid);
}

FieldReader::
~FieldReader ()
{
  finalize();
}

void FieldReader::
set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
  EKAT_REQUIRE_MSG (atm_logger, "Error! Invalid logger pointer.\n");
  m_atm_logger = atm_logger;
}

void FieldReader::
set_fields (const std::vector<Field>& fields)
{
  EKAT_REQUIRE_MSG (m_io_grid,
      "[FieldReader] Error! Grid must be set BEFORE the fields.\n");

  if (m_fields_names.size()>0) {
    EKAT_REQUIRE_MSG (m_fields_names.size()==fields.size(),
        "[FieldReader] Error! Resetting fields with a different set of fields.\n");
    for (const auto& f_new : fields) {
      EKAT_REQUIRE_MSG (m_fields_from_user->has_field(f_new.name()),
        "[FieldReader] Error! Resetting fields with a different set of fields.\n"
        " Adding extra field '" + f_new.name() + "'.\n");

      const auto& f_old = m_fields_from_user->get_field(f_new.name());
      const auto& fl_new = f_new.get_header().get_identifier().get_layout();
      const auto& fl_old = f_old.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG (fl_new==fl_old,
        "[FieldReader] Error! Resetting fields with a different set of fields.\n"
        " Changing layout of field '" + f_new.name() + "'.\n"
        "  - old layout: " + fl_old.to_string() + "\n"
        "  - new layout: " + fl_new.to_string() + "\n");
    }
  }

  m_fields_from_user = std::make_shared<FieldManager>(m_io_grid,RepoState::Closed);
  m_fields_for_scorpio = std::make_shared<FieldManager>(m_io_grid,RepoState::Closed);
  m_fields_names.clear();
  for (const auto& f : fields) {
    m_fields_names.push_back(f.name());
    m_fields_from_user->add_field(f);

    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();

    // If we can alias the field, do it,otherwise, create a clone.
    bool can_alias = fh.get_parent()==nullptr && fap.get_padding()==0;
    if (can_alias) {
      m_fields_for_scorpio->add_field(f);
    } else {
      // We have padding, or the field is a subfield (or both).
      // Either way, we need a temporary.
      // NOTE: Field::clone honors the max pack size for the allocation,
      //       which would defy one of the reason for cloning it...
      Field copy (f.get_header().get_identifier());
      copy.allocate_view();
      m_fields_for_scorpio->add_field(copy);
    }
  }

  m_fields_inited = true;
}

void FieldReader::
reset_filename (const std::string& filename,
                const std::string& iotype)
{
  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }
  m_filename = filename;
  init_scorpio_structures(iotype);
}

/* ---------------------------------------------------------- */
void FieldReader::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  EKAT_REQUIRE_MSG (grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");
  EKAT_REQUIRE_MSG (
      (grid->get_global_max_dof_gid()-grid->get_global_min_dof_gid()+1)==grid->get_num_global_dofs(),
      "Error! IO requires DOF gids to (globally)  be in interval [gid_0,gid_0+num_global_dofs).\n"
      "   - global min GID : " + std::to_string(grid->get_global_min_dof_gid()) + "\n"
      "   - global max GID : " + std::to_string(grid->get_global_max_dof_gid()) + "\n"
      "   - num global dofs: " + std::to_string(grid->get_num_global_dofs()) + "\n");

  // The grid is good. Store it.
  m_io_grid = grid;
}

/* ---------------------------------------------------------- */
// Note: The (zero-based) time_index argument provides a way to control which
//       time step to read input from in the file.  If a negative number is
//       provided the routine will read input at the last time index present in the file
void FieldReader::read_variables (const int time_index)
{
  m_atm_logger->info("[EAMxx::scorpio_input] Reading variables from file");
  m_atm_logger->info("  file name: " + m_filename);
  m_atm_logger->info("  var names: " + ekat::join(m_fields_names,", "));
  if (time_index!=-1) {
    m_atm_logger->info("  time idx : " + std::to_string(time_index));
  }

  EKAT_REQUIRE_MSG (m_fields_inited and m_scorpio_inited,
      "Error! Internal structures not fully inited yet. Did you forget to call 'init(..)'?\n");

  for (auto const& name : m_fields_names) {

    auto f_scorpio = m_fields_for_scorpio->get_field(name);
    auto f_user    = m_fields_from_user->get_field(name);

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
}

/* ---------------------------------------------------------- */
void FieldReader::finalize()
{
  if (m_scorpio_inited) {
    scorpio::release_file(m_filename);
  }

  m_fields_from_user   = nullptr;
  m_fields_for_scorpio = nullptr;
  m_io_grid            = nullptr;
  m_fields_names.clear();

  m_fields_inited  = false;
  m_scorpio_inited = false;
}

void FieldReader::init_scorpio_structures(const std::string& iotype)
{
  EKAT_REQUIRE_MSG (m_fields_inited,
      "Error! Cannot init scorpio structures until fields/views have been set.\n");

  scorpio::register_file(m_filename,scorpio::Read,scorpio::str2iotype(iotype));

  // Some input files have the "time" dimension as non-unlimited. This messes up our
  // scorpio interface, which stores a pointer to a "time" dim to be used to read/write
  // slices. This ptr is automatically inited to the unlimited dim in the file. If there is
  // no unlim dim, this ptr remains inited.
  if (not scorpio::has_time_dim(m_filename) and scorpio::has_dim(m_filename,"time")) {
    scorpio::mark_dim_as_time(m_filename,"time");
  }

  // Check variables are in the input file
  for (const auto & [name, f] : m_fields_for_scorpio->get_repo()) {
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
FieldReader::get_vec_of_dims(const FieldLayout& layout)
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
void FieldReader::set_decompositions()
{
  using namespace ShortFieldTagsNames;

  // First, check if any of the vars is indeed partitioned
  const auto decomp_tag  = m_io_grid->get_partitioned_dim_tag();

  bool has_decomposed_layouts = false;
  for (const auto & [name, f] : m_fields_for_scorpio->get_repo()) {
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
