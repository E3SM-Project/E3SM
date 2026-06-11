#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_zip.hpp>

#include <memory>
#include <numeric>

namespace scream
{

FieldReader::
~FieldReader ()
{
  clean_up ();
}

void FieldReader::
set_fields (const std::vector<Field>& fields)
{
  m_fields = fields;
  m_reader_state |= NEW_FIELDS;
}

void FieldReader::
set_file_specs (const std::string& filename,
                const std::string& io_type)
{
  set_file_specs(filename,{},io_type);
}

void FieldReader::
set_file_specs (const std::string& filename,
                const std::map<std::string,std::string>& tag_rename,
                const std::string& io_type)

{
  EKAT_REQUIRE_MSG (filename!="",
      "[FieldReader::set_file_specs] Error! Input filename is empty.\n");

  if (filename!=m_filename) {
    if (m_filename!="") {
      scorpio::release_file(m_filename);
    }

    m_filename = filename;
    m_tag_rename = tag_rename;

    scorpio::register_file(m_filename,scorpio::FileMode::Read,scorpio::str2iotype(io_type));

    // Some input files have the "time" dimension as non-unlimited. This can sometime mess up
    // our interfaces, which stores a pointer to a "time" dim to be used to read/write
    // the desired time slices. This ptr is automatically inited to the unlimited dim in the file.
    // If there is no unlimited dim, this ptr remains uninited, which means that the
    // eamxx scorpio interface will treat "time" like any other dim (yielding n+1-dim fields)
    if (not scorpio::has_time_dim(m_filename) and scorpio::has_dim(m_filename,"time")) {
      scorpio::mark_dim_as_time(m_filename,"time");
    }

    m_reader_state |= NEW_FILE;
  } else {
    EKAT_REQUIRE_MSG (tag_rename==m_tag_rename,
        "[FieldReader::set_file_specs] Error! Filename did not change, but tag_rename did.\n"
        " - file name: " + filename + "\n");
  }
}

void FieldReader::read (const int time_index)
{
  if (m_reader_state!=CLEAN)
    setup_internals ();

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
  m_tag_rename = {};

  m_dim_decomp_offsets = {};
  m_dim_decomp_name = "";

  m_reader_state = CLEAN;
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
  m_dim_decomp_offsets.resize(local_dim);
  for (int idof=0; idof<local_dim; ++idof) {
    m_dim_decomp_offsets[idof] = gids_h[idof] - min_gid;
  }
  m_dim_decomp_name = layout.name(0);

  m_reader_state |= NEW_DECOMP;
}

void FieldReader::setup_internals ()
{
  if (m_dim_decomp_name!="" and m_reader_state & (NEW_FILE | NEW_DECOMP)) {
    if (std::is_same_v<std::int64_t,scorpio::offset_t>) {
      scorpio::set_dim_decomp(m_filename, m_dim_decomp_name, m_dim_decomp_offsets);
    } else {
      std::vector<scorpio::offset_t> offsets(m_dim_decomp_offsets.begin(),
                                             m_dim_decomp_offsets.end());
      scorpio::set_dim_decomp(m_filename, m_dim_decomp_name, offsets);
    }
  }

  if (m_reader_state & (NEW_FIELDS | NEW_FILE)) {
    // First, check fields are in the file, with correct layout
    for (const auto & f : m_fields) {
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
        if (m_tag_rename.count(io_dims[i])>0) {
          io_dims[i] = m_tag_rename.at(io_dims[i]);
        }
      }

      const auto& var = scorpio::get_var(m_filename,f.name());
      EKAT_REQUIRE_MSG (var.dim_names()==io_dims,
          "Error! Layout mismatch for input file variable.\n"
          " - filename: " + m_filename + "\n"
          " - varname : " + f.name() + "\n"
          " - expected dim names: " + ekat::join(io_dims,",") + "\n"
          " - dims from file    : " + ekat::join(var.dim_names(),",") + "\n");

      for (int i=0; i<layout.rank(); ++i) {
        if (io_dims[i]!=m_dim_decomp_name) {
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

    // Now create io_fields
    // IO fields MUST be contiguous, so if padded or with parent, create a new one.
    // Otherwise, we can alias input field, possibly with renamed tags
    m_io_fields = {};
    auto io_map = m_layout_to_io_field;
    m_layout_to_io_field.clear();
    for (const auto& f : m_fields) {
      const auto& fh  = f.get_header();
      const auto& fap = fh.get_alloc_properties();

      const auto& fid    = fh.get_identifier();
      const auto& layout = fid.get_layout();
      auto key = e2str(fid.data_type()) + "_" + ekat::join(layout.dims(),"_");
      if (fh.get_parent() or fap.get_padding()>0) {
        auto it = io_map.try_emplace(key,fid.clone(key),true);
        m_io_fields.push_back(it.first->second.alias(f.name(),m_tag_rename));
        m_layout_to_io_field[key] = it.first->second;
      } else {
        m_io_fields.push_back(f.alias(f.name(),m_tag_rename));
        m_layout_to_io_field[key] = m_io_fields.back();
      }
    }
  }

  m_reader_state = CLEAN;
}

// ------ Free functions ------- //

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const int time_index)
{
  // FieldReader reader(filename,fields);
  FieldReader reader;
  reader.set_file_specs(filename);
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
  FieldReader reader;
  reader.set_file_specs(filename);
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
