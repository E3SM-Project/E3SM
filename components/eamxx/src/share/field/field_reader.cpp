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
  // IO fields MUST be contiguous, so if padded or with parent, create a new one
  auto make_io_field = [](const Field& f) {
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    if (fh.get_parent() or fap.get_padding()>0) {
      auto fid = fh.get_identifier().clone();
      return Field(fid,true);
    }

    return f;
  };

  bool fields_changed = true;
  if (m_fields.size()>0) {
    // If resetting the fields, ensure we are not changing the field,
    // as the scorpio structures have already been inited
    EKAT_REQUIRE_MSG (fields.size()==m_fields.size(),
        "[FieldReader::set_fields] Error! Cannot change the number of fields once set.\n");

    // Assume fields did not change. If 1+ did, we'll flip this to true
    fields_changed = false;

    auto find_f = [&](const auto& name, auto& fs) {
      std::pair<bool,Field*> p = {false,nullptr};
      for (auto& f : fs)
        if (f.name()==name) {
          p.first = true;
          p.second = &f;
        }
      return p;
    };
    int pos = 0;
    for (const auto& f_new : fields) {
      auto [found,f_old] = find_f(f_new.name(),m_fields);
      EKAT_REQUIRE_MSG(found,
        "[FieldReader::set_fields] Error! Cannot change field names.\n"
        " - new field name: " + f_new.name() + "\n"
        " - stored fields : " + ekat::join(m_fields,[](auto f){return f.name();},",") + "\n");

      // Check Layouts
      const auto& l_old = f_old->get_header().get_identifier().get_layout();
      const auto& l_new = f_new.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG(l_old==l_new,
        "[FieldReader::set_fields] Error! Cannot change field layout.\n"
        " - field name: " + f_new.name() + "\n"
        " - new layout: " + l_new.to_string() + "\n"
        " - old layout: " + l_old.to_string() + "\n");

      if (not f_new.is_aliasing(*f_old)) {
        fields_changed = true;
        m_io_fields[pos] = make_io_field(f_new);
        m_fields[pos] = f_new;
      }
      ++pos;
    }
    return;
  }

  m_fields = fields;
  m_io_fields.reserve(fields.size());
  for (const auto& f : fields) {
    m_io_fields.push_back(make_io_field(f));
  }

  if (m_inited and fields_changed)
    init_scorpio_structures();
}

void FieldReader::
set_file_specs (const std::string& filename, const std::string& io_type)
{
  set_file_specs(filename,{},io_type);
}

void FieldReader::
set_file_specs (const std::string& filename, const std::map<FieldTag,std::string>& tag2name, const std::string& io_type)
{
  EKAT_REQUIRE_MSG (filename!="",
      "[FieldReader::set_file_specs] Error! Input filename is empty.\n");

  // If resetting to same file, there's nothing to do
  if (filename==m_filename and tag2name==m_tag2name)
    return;

  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }
  m_filename = filename;
  m_tag2name = tag2name;

  scorpio::register_file(filename,scorpio::Read,scorpio::str2iotype(io_type));

  // Some input files have the "time" dimension as non-unlimited. This messes up our
  // scorpio interface, which stores a pointer to a "time" dim to be used to read/write
  // slices. This ptr is automatically inited to the unlimited dim in the file. If there is
  // no unlim dim, this ptr remains inited.
  if (not scorpio::has_time_dim(m_filename) and scorpio::has_dim(m_filename,"time")) {
    scorpio::mark_dim_as_time(m_filename,"time");
  }

  if (m_inited) {
    init_scorpio_structures();
  }
}

void FieldReader::read (const int time_index)
{
  EKAT_REQUIRE_MSG (m_filename!="",
      "[FieldReader::read] Error! The class was not correctly initialized.\n"
      " Did you forget to call 'set_file_specs(..)'?\n");

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
    if (not f_io.is_aliasing(f)) {
      f.deep_copy(f_io);
    }
  }
}

void FieldReader::clean_up ()
{
  if (m_filename!="")
    scorpio::release_file(m_filename);

  m_fields = {};
  m_io_fields = {};
  m_decomp_dim_gids = {};
  m_filename = "";
}

void FieldReader::init_scorpio_structures()
{
  EKAT_REQUIRE_MSG (m_filename!="",
      "[FieldReader::init_scorpio_structures] Error! Filename not set.\n"
      " Did you forget to call 'set_filename'?\n");

  // If a decomp was set, process it
  if (m_decomp_dim_gids.is_allocated()) {
    const auto& gids = m_decomp_dim_gids;
    const auto& layout = gids.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG (gids.rank()==1,
        "[FieldReader::set_dim_decomp] Error! GIDs field must have rank 1.\n"
        " - field name: " + gids.name() + "\n"
        " - field layout: " + layout.to_string() + "\n");
    EKAT_REQUIRE_MSG (gids.data_type()==DataType::IntType,
        "[FieldReader::set_dim_decomp] Error! GIDs field must have data type IntType.\n"
        " - field name : " + gids.name() + "\n"
        " - field dtype: " + e2str(gids.data_type()) + "\n");

    auto decomp_tag = layout.tag(0);
    bool has_decomposed_layouts = false;
    for (const auto & f : m_io_fields) {
      const auto& f_layout = f.get_header().get_identifier().get_layout();
      has_decomposed_layouts |= f_layout.has_tag(decomp_tag);
    }

    if (has_decomposed_layouts) {
      // Set the decomposition for the partitioned dimension
      const int local_dim = layout.size();
      auto gids_h = gids.get_view<const int*,Host>();
      auto min_gid = field_min(gids,&m_comm).as<int>();
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
  }

  // Check variables are in the input file, with the correct dims
  auto decomp_tag = m_decomp_dim_gids.is_allocated()
                  ? m_decomp_dim_gids.get_header().get_identifier().get_layout().tag(0)
                  : FieldTag::Invalid;

  for (const auto & f : m_io_fields) {
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

    // Check that the variable is in the file.
    EKAT_REQUIRE_MSG (scorpio::has_var(m_filename,f.name()),
        "Error! Input file does not store a required variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n");

    const auto& var = scorpio::get_var(m_filename,f.name());
    EKAT_REQUIRE_MSG (var.dim_names()==io_dims,
        "Error! Layout mismatch for input file variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n"
        " - expected dims : " + ekat::join(io_dims,",") + "\n"
        " - dims from file: " + ekat::join(var.dim_names(),",") + "\n");

    for (int i=0; i<layout.rank(); ++i) {
      if (layout.tag(i)==decomp_tag)
        continue;

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

    // Ensure that we can read the var using Real data type
    scorpio::change_var_dtype (m_filename,f.name(),"real");
  }

  // From now on, calling set_file_specs triggers a re-initialization
  m_inited = true;
}

void FieldReader::set_dim_decomp(const Field& gids, const ekat::Comm& comm)
{
  m_decomp_dim_gids = gids;
  m_comm = comm;
}

// ------ Free functions ------- //

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const int time_index)
{
  // FieldReader reader(filename,fields);
  FieldReader reader;
  reader.set_fields(fields);
  reader.set_file_specs(filename);
  reader.init_scorpio_structures();
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
  reader.init_scorpio_structures();
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
