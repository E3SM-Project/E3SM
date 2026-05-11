#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_zip.hpp>

#include <memory>
#include <numeric>

namespace scream
{

namespace {
std::string get_io_dtype (const Field& f)
{
  std::string io_dt = "invalid";
  switch (f.data_type()) {
    case DataType::IntType:
      io_dt = "int";
      break;
    case DataType::FloatType:
      io_dt = "float";
      break;
    case DataType::DoubleType:
      io_dt = "double";
      break;
    default:
      EKAT_ERROR_MSG (
          "Error! Unsupported/unrecognized data type.\n"
          " - field name: " + f.name() + "\n");
  }
  return io_dt;
}
} // anonymous namespace

FieldReader::
FieldReader (const std::string& filename, const std::string& io_type)
 : FieldReader(filename,{},io_type)
{
  // Noting to do here
}

FieldReader::
FieldReader (const std::string& filename,
             const std::map<FieldTag,std::string>& tag2name,
             const std::string& io_type)
{
  set_file_specs(filename,tag2name,io_type);
}

FieldReader::
~FieldReader ()
{
  clean_up ();
}

void FieldReader::
set_fields (const std::vector<Field>& fields)
{
  check_fields(fields);

  m_fields = fields;
  m_io_fields.clear();
  m_io_fields.reserve(fields.size());
  for (const auto& f : fields) {
    // Store the io field
    m_io_fields.push_back(make_io_field(f));
  }
}

void FieldReader::
set_file_specs (const std::string& filename,
                const std::string& io_type)
{
  set_file_specs(filename,{},io_type);
}

void FieldReader::
set_file_specs (const std::string& filename,
                const std::map<FieldTag,std::string>& tag2name,
                const std::string& io_type)

{
  EKAT_REQUIRE_MSG (filename!="",
      "[FieldReader::set_file_specs] Error! Input filename is empty.\n");

  if (filename==m_filename) {
    // Ensure that the vars dtype is correct (in case dtype!=nc_dtype
    for (const auto& f : m_io_fields) {
      auto dtype = get_io_dtype(f);
      scorpio::change_var_dtype(m_filename,f.name(),dtype);
    }
    return;
  }

  m_tag2name = tag2name;

  scorpio::register_file(filename,scorpio::FileMode::Read,scorpio::str2iotype(io_type));

  // Some input files have the "time" dimension as non-unlimited. This can sometime mess up
  // our interfaces, which stores a pointer to a "time" dim to be used to read/write
  // the desired time slices. This ptr is automatically inited to the unlimited dim in the file.
  // If there is no unlimited dim, this ptr remains uninited, which means that the
  // eamxx scorpio interface will treat "time" like any other dim (yielding n+1-dim fields)
  if (not scorpio::has_time_dim(filename) and scorpio::has_dim(filename,"time")) {
    scorpio::mark_dim_as_time(filename,"time");
  }

  if (m_filename!="") {
    // Transfer dim decompositions
    auto dim_decomps = scorpio::get_dim_decomps(m_filename);
    for (auto [name,decomp] : dim_decomps) {
      if (decomp->dims.size()==1) {
        scorpio::set_dim_decomp(filename,decomp->dims[0]->name,decomp->offsets);
      } else {
        std::vector<std::string> dimnames;
        for (auto d : decomp->dims)
          dimnames.push_back(d->name);
        scorpio::set_dims_decomp(filename,dimnames,decomp->offsets);
      }
    }
  }

  if (m_fields.size()>0) {
    // If the new file stores fields with same dtype, make_io_field would return the field itself
    check_fields(m_fields);
    for (auto& f : m_io_fields) {
      f = make_io_field(f);
      auto dtype = get_io_dtype(f);
      scorpio::change_var_dtype(m_filename,f.name(),dtype);
    }
  }

  if (m_filename!="") {
    // Release old file (if any) and update stored file name
    scorpio::release_file(m_filename);
  }
  m_filename = filename;
}

void FieldReader::read (const int time_index)
{
  for (const auto& [f, f_io] : ekat::zip(m_fields,m_io_fields)) {
    switch (f_io.data_type()) {
      case DataType::DoubleType:
        scorpio::read_var(m_filename,f.name(),f_io.get_internal_view_data<double,Host>(),time_index);
        break;
      case DataType::FloatType:
        scorpio::read_var(m_filename,f.name(),f_io.get_internal_view_data<float,Host>(),time_index);
        break;
      case DataType::IntType:
        scorpio::read_var(m_filename,f.name(),f_io.get_internal_view_data<int,Host>(),time_index);
        break;
      default:
        EKAT_ERROR_MSG (
            "Error! Unsupported/unrecognized data type while reading field from file.\n"
            " - file name : " + m_filename + "\n"
            " - field name: " + f.name() + "\n");
    }

    f_io.sync_to_dev();
    f.deep_copy(f_io); // If f_io is aliasing f, this is a no-op
  }
}

void FieldReader::clean_up ()
{
  if (m_filename!="")
    scorpio::release_file(m_filename);

  m_fields = {};
  m_io_fields = {};
  m_filename = "";
}

void FieldReader::set_dim_decomp(const Field& gids, const ekat::Comm& comm)
{
  const auto& layout = gids.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (gids.is_allocated(),
      "[FieldReader::set_dim_decomp] Error! Invalid gids field.\n"
      " - gids field name: " + gids.name() + "\n");
  EKAT_REQUIRE_MSG (gids.rank()==1,
      "[FieldReader::set_dim_decomp] Error! GIDs field must have rank 1.\n"
      " - field name: " + gids.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG (gids.data_type()==DataType::IntType,
      "[FieldReader::set_dim_decomp] Error! GIDs field must have data type IntType.\n"
      " - field name : " + gids.name() + "\n"
      " - field dtype: " + e2str(gids.data_type()) + "\n");

  // Set the decomposition for the partitioned dimension
  const int local_dim = layout.size();
  auto gids_h = gids.get_view<const int*,Host>();
  auto min_gid = field_min(gids,&comm).as<int>();
  std::vector<scorpio::offset_t> offsets(local_dim);
  for (int idof=0; idof<local_dim; ++idof) {
    offsets[idof] = gids_h[idof] - min_gid;
  }
  if (m_tag2name.count(layout.tag(0))>0) {
    scorpio::set_dim_decomp(m_filename,m_tag2name.at(layout.tag(0)),offsets);
  } else {
    scorpio::set_dim_decomp(m_filename,layout.name(0),offsets);
  }
}

Field FieldReader::make_io_field (const Field& f)
{
  // IO fields MUST be contiguous, so if padded or with parent, create a new one.
  // Also, if the data type is different, create a clone to match the file data type
  // NOTE: we can use the same io field for multiple input fields, if dtype and layout are the same
  const auto& fh  = f.get_header();
  const auto& fap = fh.get_alloc_properties();

  if (fh.get_parent() or fap.get_padding()>0) {
    const auto& fid = fh.get_identifier();
    auto it = m_layout_to_io_field.try_emplace(fid.get_layout().to_string(),f.clone());
    return it.first->second;
  } else {
    return f.alias(f.name(),m_tag2name);
  }
}

void FieldReader::check_fields(const std::vector<Field>& fields) const
{
  // Check variables are in the input file, with the correct dims
  for (const auto & f : fields) {
    // Check that the variable is in the file.
    EKAT_REQUIRE_MSG (scorpio::has_var(m_filename,f.name()),
        "Error! Input file does not store a required variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n");

    const auto& layout = f.get_header().get_identifier().get_layout();
    auto io_dims = layout.names();
    for (int i=0; i<layout.rank(); ++i) {
      if (io_dims[i]=="dim")
        io_dims[i] += std::to_string(layout.dim(i));
      auto t = layout.tag(i);
      if (m_tag2name.count(t)>0) {
        io_dims[i] = m_tag2name.at(t);
      }
    }

    const auto& var = scorpio::get_var(m_filename,f.name());
    EKAT_REQUIRE_MSG (var.dim_names()==io_dims,
        "Error! Layout mismatch for input file variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n"
        " - expected dims : " + ekat::join(io_dims,",") + "\n"
        " - dims from file: " + ekat::join(var.dim_names(),",") + "\n");

    for (int i=0; i<layout.rank(); ++i) {
      const int dim_len  = scorpio::get_dimlen_local(m_filename,io_dims[i]);
      const int dim_glen  = scorpio::get_dimlen(m_filename,io_dims[i]);
      if (dim_len==dim_glen) {
        const int file_len  = scorpio::get_dimlen(m_filename,io_dims[i]);
        EKAT_REQUIRE_MSG (layout.dim(i)==file_len,
            "Error! Dimension mismatch for input file variable.\n"
          " - filename : " + m_filename + "\n"
          " - varname  : " + f.name() + "\n"
          " - var dims : " + ekat::join(io_dims,",") + "\n"
          " - dim name : " + io_dims[i] + "\n"
          " - expected extent : " + std::to_string(layout.dim(i)) + "\n"
          " - extent from file: " + std::to_string(file_len) + "\n");
      }
    }
  }
}

// ------ Free functions ------- //

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const int time_index)
{
  // FieldReader reader(filename,fields);
  FieldReader reader(filename);
  reader.set_fields(fields);
  reader.read(time_index);
}
void read_fields (const std::string& filename,
                  const std::initializer_list<Field>& fields,
                  const int time_index)
{
  read_fields(filename,std::vector<Field>(fields),time_index);
}
void read_fields (const std::string& filename,
                  const std::map<std::string,Field>& fields,
                  const int time_index)
{
  std::vector<Field> fields_vec;
  for (const auto& [name,f] : fields)
    fields_vec.push_back(f.alias(name));
  read_fields(filename,fields_vec,time_index);
}

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index)
{
  // FieldReader reader(filename,fields,decomp_gids,comm);
  FieldReader reader(filename);
  reader.set_dim_decomp(decomp_gids,comm);
  reader.set_fields(fields);
  reader.read(time_index);
}
void read_fields (const std::string& filename,
                  const std::initializer_list<Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index)
{
  read_fields(filename,std::vector<Field>(fields),decomp_gids,comm,time_index);
}
void read_fields (const std::string& filename,
                  const std::map<std::string,Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index)
{
  std::vector<Field> fields_vec;
  for (const auto& [name,f] : fields)
    fields_vec.push_back(f.alias(name));
  read_fields(filename,fields_vec,decomp_gids,comm,time_index);
}

} // namespace scream
