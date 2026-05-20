#include "share/field/field_writer.hpp"

#include "share/field/field_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_zip.hpp>

#include <filesystem>
#include <numeric>

namespace scream
{

FieldWriter::
~FieldWriter ()
{
  clean_up();
}

void FieldWriter::
set_fields (const std::vector<Field>& fields)
{
  auto make_io_field = [](const Field& f) {
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    if (fh.get_parent() or fap.get_padding()>0) {
      auto fid = fh.get_identifier().clone();
      return Field(fid,true);
    }

    return f;
  };

  if (m_fields.size()>0) {
    EKAT_REQUIRE_MSG (fields.size()==m_fields.size(),
        "[FieldWriter::set_fields] Error! Cannot change the number of fields once set.\n");

    auto find_f = [&](const auto& name, auto& fs) {
      std::pair<bool,Field*> p = {false,nullptr};
      for (auto& f : fs) {
        if (f.name()==name) {
          p.first = true;
          p.second = &f;
        }
      }
      return p;
    };

    int pos = 0;
    for (const auto& f_new : fields) {
      auto [found,f_old] = find_f(f_new.name(),m_fields);
      EKAT_REQUIRE_MSG(found,
          "[FieldWriter::set_fields] Error! Cannot change field names.\n"
          " - new field name: " + f_new.name() + "\n"
          " - stored fields : " + ekat::join(m_fields,[](auto f){return f.name();},",") + "\n");

      const auto& l_old = f_old->get_header().get_identifier().get_layout();
      const auto& l_new = f_new.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG(l_old==l_new,
          "[FieldWriter::set_fields] Error! Cannot change field layout.\n"
          " - field name: " + f_new.name() + "\n"
          " - new layout: " + l_new.to_string() + "\n"
          " - old layout: " + l_old.to_string() + "\n");

      if (not f_new.is_aliasing(*f_old)) {
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

  if (m_inited) {
    init_scorpio_structures();
  }
}

void FieldWriter::
set_file_specs (const std::string& filename, const std::string& io_type)
{
  set_file_specs(filename,{},io_type);
}

void FieldWriter::
set_file_specs (const std::string& filename,
                const std::map<FieldTag,std::string>& tag2name,
                const std::string& io_type)
{
  EKAT_REQUIRE_MSG (filename!="",
      "[FieldWriter::set_file_specs] Error! Output filename is empty.\n");

  if (filename==m_filename and tag2name==m_tag2name) {
    return;
  }

  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }

  m_filename = filename;
  m_tag2name = tag2name;

  const auto exists = std::filesystem::exists(filename);
  const auto mode = exists ? scorpio::Append : scorpio::Write;
  scorpio::register_file(filename,mode,scorpio::str2iotype(io_type));

  if (m_inited) {
    init_scorpio_structures();
  }
}

void FieldWriter::
set_dim_decomp(const Field& gids, const ekat::Comm& comm)
{
  m_decomp_dim_gids = gids;
  m_comm = comm;
}

void FieldWriter::
set_time_dependent (const bool time_dep,
                    const std::string& time_units,
                    const std::string& time_name)
{
  EKAT_REQUIRE_MSG (time_name!="",
      "[FieldWriter::set_time_dependent] Error! Invalid empty time dimension name.\n");

  if (m_inited) {
    EKAT_REQUIRE_MSG (m_time_dep==time_dep,
        "[FieldWriter::set_time_dependent] Error! Cannot change time dependency after initialization.\n");
    if (m_time_dep) {
      EKAT_REQUIRE_MSG (m_time_name==time_name,
          "[FieldWriter::set_time_dependent] Error! Cannot change time dimension name after initialization.\n");
    }
    return;
  }

  m_time_dep = time_dep;
  m_time_units = time_units;
  m_time_name = time_name;
}

void FieldWriter::
init_scorpio_structures()
{
  EKAT_REQUIRE_MSG (m_filename!="",
      "[FieldWriter::init_scorpio_structures] Error! Filename not set.\n"
      " Did you forget to call 'set_file_specs' ?\n");

  EKAT_REQUIRE_MSG (m_fields.size()>0,
      "[FieldWriter::init_scorpio_structures] Error! No fields registered.\n"
      " Did you forget to call 'set_fields' ?\n");

  bool in_def_phase = false;
  auto ensure_def_phase = [&]() {
    if (not in_def_phase) {
      scorpio::redef(m_filename);
      in_def_phase = true;
    }
  };

  auto decomp_tag = m_decomp_dim_gids.is_allocated()
                  ? m_decomp_dim_gids.get_header().get_identifier().get_layout().tag(0)
                  : FieldTag::Invalid;

  if (m_decomp_dim_gids.is_allocated()) {
    const auto& gids = m_decomp_dim_gids;
    const auto& layout = gids.get_header().get_identifier().get_layout();

    EKAT_REQUIRE_MSG (gids.rank()==1,
        "[FieldWriter::init_scorpio_structures] Error! GIDs field must have rank 1.\n"
        " - field name  : " + gids.name() + "\n"
        " - field layout: " + layout.to_string() + "\n");
    EKAT_REQUIRE_MSG (gids.data_type()==DataType::IntType,
        "[FieldWriter::init_scorpio_structures] Error! GIDs field must have data type IntType.\n"
        " - field name : " + gids.name() + "\n"
        " - field dtype: " + e2str(gids.data_type()) + "\n");
  }

  auto get_dim_name = [&](const FieldLayout& layout, const int i) {
    auto dim_name = layout.name(i);
    if (dim_name=="dim") {
      dim_name += std::to_string(layout.dim(i));
    }

    const auto t = layout.tag(i);
    if (m_tag2name.count(t)>0) {
      dim_name = m_tag2name.at(t);
    }
    return dim_name;
  };

  for (const auto& f : m_io_fields) {
    const auto& layout = f.get_header().get_identifier().get_layout();
    for (int i=0; i<layout.rank(); ++i) {
      const auto dim_name = get_dim_name(layout,i);

      int dim_len = layout.dim(i);
      if (layout.tag(i)==decomp_tag) {
        EKAT_REQUIRE_MSG (m_decomp_dim_gids.is_allocated(),
            "[FieldWriter::init_scorpio_structures] Error! Decomposed field found but no decomposition gids were set.\n"
            " - field name: " + f.name() + "\n"
            " - dim name  : " + dim_name + "\n");

        auto min_gid = field_min(m_decomp_dim_gids,&m_comm).as<int>();
        auto max_gid = field_max(m_decomp_dim_gids,&m_comm).as<int>();
        dim_len = max_gid-min_gid+1;
      }

      if (scorpio::has_dim(m_filename,dim_name)) {
        if (layout.tag(i)!=decomp_tag) {
          const int file_len = scorpio::get_dimlen(m_filename,dim_name);
          EKAT_REQUIRE_MSG (file_len==dim_len,
              "[FieldWriter::init_scorpio_structures] Error! Dimension mismatch in output file.\n"
              " - filename: " + m_filename + "\n"
              " - field   : " + f.name() + "\n"
              " - dim name: " + dim_name + "\n"
              " - expected len: " + std::to_string(dim_len) + "\n"
              " - file len    : " + std::to_string(file_len) + "\n");
        }
      } else {
        ensure_def_phase();
        scorpio::define_dim(m_filename,dim_name,dim_len);
      }
    }
  }

  if (m_time_dep and not scorpio::has_time_dim(m_filename)) {
    ensure_def_phase();
    scorpio::define_time(m_filename,m_time_units,m_time_name);
  }

  if (m_decomp_dim_gids.is_allocated()) {
    const auto& gids = m_decomp_dim_gids;
    const auto& layout = gids.get_header().get_identifier().get_layout();

    bool has_decomposed_layouts = false;
    for (const auto & f : m_io_fields) {
      const auto& f_layout = f.get_header().get_identifier().get_layout();
      has_decomposed_layouts |= f_layout.has_tag(layout.tag(0));
    }

    if (has_decomposed_layouts) {
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

  for (const auto& f : m_io_fields) {
    const auto& fid = f.get_header().get_identifier();
    const auto& layout = fid.get_layout();

    std::vector<std::string> io_dims(layout.rank());
    for (int i=0; i<layout.rank(); ++i) {
      io_dims[i] = get_dim_name(layout,i);
    }

    std::string dtype;
    switch (f.data_type()) {
      case DataType::IntType:
        dtype = "int";
        break;
      case DataType::FloatType:
        dtype = "float";
        break;
      case DataType::DoubleType:
        dtype = "double";
        break;
      default:
        EKAT_ERROR_MSG (
            "Error! Unsupported/unrecognized data type while defining output field.\n"
            " - file name : " + m_filename + "\n"
            " - field name: " + f.name() + "\n");
    }

    ensure_def_phase();
    scorpio::define_var(m_filename,f.name(),fid.get_units().to_string(),io_dims,dtype,dtype,m_time_dep);
  }

  if (in_def_phase) {
    scorpio::enddef(m_filename);
  }

  m_inited = true;
}

void FieldWriter::
write_impl ()
{
  EKAT_REQUIRE_MSG (m_filename!="",
      "[FieldWriter::write] Error! The class was not correctly initialized.\n"
      " Did you forget to call 'set_file_specs(..)'?\n");

  for (const auto& [f,f_io] : ekat::zip(m_fields,m_io_fields)) {
    if (not f_io.is_aliasing(f)) {
      f_io.deep_copy(f);
    }
    f_io.sync_to_host();

    switch (f_io.data_type()) {
      case DataType::DoubleType:
        scorpio::write_var(m_filename,f.name(),f_io.get_internal_view_data<const double,Host>());
        break;
      case DataType::FloatType:
        scorpio::write_var(m_filename,f.name(),f_io.get_internal_view_data<const float,Host>());
        break;
      case DataType::IntType:
        scorpio::write_var(m_filename,f.name(),f_io.get_internal_view_data<const int,Host>());
        break;
      default:
        EKAT_ERROR_MSG (
            "Error! Unsupported/unrecognized data type while writing field to file.\n"
            " - file name : " + m_filename + "\n"
            " - field name: " + f.name() + "\n");
    }
  }
}

void FieldWriter::
write ()
{
  EKAT_REQUIRE_MSG (not m_time_dep,
      "[FieldWriter::write] Error! Time-dependent writer requires write(time).\n");

  write_impl();
}

void FieldWriter::
write (const double time)
{
  EKAT_REQUIRE_MSG (m_time_dep,
      "[FieldWriter::write(time)] Error! Writer was initialized as time-independent.\n");

  scorpio::update_time(m_filename,time);
  write_impl();
}

void FieldWriter::
clean_up ()
{
  if (m_filename!="") {
    scorpio::release_file(m_filename);
  }

  m_fields = {};
  m_io_fields = {};
  m_decomp_dim_gids = {};
  m_filename = "";
  m_inited = false;
}

void write_fields (const std::string& filename,
                   const std::vector<Field>& fields)
{
  FieldWriter writer;
  writer.set_fields(fields);
  writer.set_file_specs(filename);
  writer.init_scorpio_structures();
  writer.write();
}

void write_fields (const std::string& filename,
                   const std::initializer_list<Field>& fields)
{
  write_fields(filename,std::vector<Field>(fields));
}

void write_fields (const std::string& filename,
                   const std::map<std::string,Field>& fields)
{
  std::vector<Field> fields_vec;
  for (const auto& [name,f] : fields) {
    fields_vec.push_back(f.alias(name));
  }
  write_fields(filename,fields_vec);
}

void write_fields (const std::string& filename,
                   const std::vector<Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm)
{
  FieldWriter writer;
  writer.set_file_specs(filename);
  writer.set_dim_decomp(decomp_gids,comm);
  writer.set_fields(fields);
  writer.init_scorpio_structures();
  writer.write();
}

void write_fields (const std::string& filename,
                   const std::initializer_list<Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm)
{
  write_fields(filename,std::vector<Field>(fields),decomp_gids,comm);
}

void write_fields (const std::string& filename,
                   const std::map<std::string,Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm)
{
  std::vector<Field> fields_vec;
  for (const auto& [name,f] : fields) {
    fields_vec.push_back(f.alias(name));
  }
  write_fields(filename,fields_vec,decomp_gids,comm);
}

} // namespace scream
