#include "share/field/field_utils.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_assert.hpp>
#include <ekat_std_utils.hpp>

#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using namespace scream;

using decomp_key_t = std::pair<std::string,FieldTag>;
using stratts_t = std::map<std::string,std::string>;

struct DecompInfo {
  std::string dim_name;
  int         global_extent;
  std::vector<scorpio::offset_t> offsets;
};

// std::vector<Field>
// group_to_fields (const FieldGroup& group)
// {
//   std::vector<Field> fields;
//   fields.reserve(group.m_individual_fields.size());

//   for (const auto& it : group.m_individual_fields) {
//     const auto& f = *it.second;
//     const auto& p = f.get_header().get_parent();
//     if (p and group.m_individual_fields.count(p->get_identifier().name())>0) {
//       continue;
//     }
//     fields.push_back(f);
//   }

//   return fields;
// }

// std::string
// get_dtype_string (const DataType dt)
// {
//   switch (dt) {
//     case DataType::IntType:
//       return "int";
//     case DataType::FloatType:
//       return "float";
//     case DataType::DoubleType:
//       return "double";
//     case DataType::RealType:
//       return "real";
//     default:
//       EKAT_ERROR_MSG ("Error! Unsupported data type in write_fields.\n");
//   }
//   return "";
// }

void
write_field_data (const std::string& filename,
                  const Field& field)
{
  Field io_field;
  const auto& fap = field.get_header().get_alloc_properties();
  if (field.get_header().get_parent()==nullptr and fap.get_padding()==0) {
    io_field = field;
  } else {
    io_field = Field(field.get_header().get_identifier());
    io_field.allocate_view();
    io_field.deep_copy(field);
  }

  io_field.sync_to_host();
  switch (io_field.data_type()) {
    case DataType::IntType:
      scorpio::write_var(filename,io_field.name(),io_field.get_internal_view_data<int,Host>());
      break;
    case DataType::FloatType:
      scorpio::write_var(filename,io_field.name(),io_field.get_internal_view_data<float,Host>());
      break;
    case DataType::DoubleType:
      scorpio::write_var(filename,io_field.name(),io_field.get_internal_view_data<double,Host>());
      break;
    case DataType::RealType:
      scorpio::write_var(filename,io_field.name(),io_field.get_internal_view_data<Real,Host>());
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported data type in write_fields.\n");
  }
}

// std::map<decomp_key_t,DecompInfo>
// create_decomp_map (const ekat::Comm* comm,
//                    const std::vector<Field>& gids_fields)
// {
//   std::map<decomp_key_t,DecompInfo> info;
//   if (comm==nullptr) {
//     return info;
//   }

//   using gid_type = AbstractGrid::gid_type;
//   for (const auto& gids : gids_fields) {
//     const auto& id = gids.get_header().get_identifier();
//     const auto& layout = id.get_layout();

//     EKAT_REQUIRE_MSG (layout.rank()==1,
//         "Error! Gids fields passed to write_fields must be rank-1.\n"
//         " - field name: " + gids.name() + "\n"
//         " - layout    : " + layout.to_string() + "\n");
//     EKAT_REQUIRE_MSG (gids.data_type()==DataType::IntType,
//         "Error! Gids fields passed to write_fields must store integer gids.\n"
//         " - field name: " + gids.name() + "\n");

//     const auto tag = layout.tag(0);
//     const auto key = std::make_pair(id.get_grid_name(),tag);
//     EKAT_REQUIRE_MSG (info.count(key)==0,
//         "Error! Duplicate gids field passed to write_fields for the same grid/tag.\n"
//         " - grid: " + id.get_grid_name() + "\n"
//         " - tag : " + e2str(tag) + "\n");

//     gids.sync_to_host();
//     auto gids_h = gids.get_view<const gid_type*,Host>();
//     const auto local_count = layout.dim(0);

//     gid_type local_min = std::numeric_limits<gid_type>::max();
//     gid_type local_max = std::numeric_limits<gid_type>::lowest();
//     for (int i=0; i<local_count; ++i) {
//       local_min = std::min(local_min,gids_h(i));
//       local_max = std::max(local_max,gids_h(i));
//     }

//     gid_type global_min = local_min;
//     gid_type global_max = local_max;
//     comm->all_reduce(&local_min,&global_min,1,MPI_MIN);
//     comm->all_reduce(&local_max,&global_max,1,MPI_MAX);
//     int global_count = 0;
//     comm->all_reduce(&local_count,&global_count,1,MPI_SUM);

//     if (global_count==0) {
//       global_min = 0;
//       global_max = -1;
//     }

//     DecompInfo decomp;
//     decomp.dim_name = layout.name(0);
//     decomp.global_extent = global_max>=global_min ? global_max-global_min+1 : 0;
//     decomp.offsets.resize(local_count);
//     for (int i=0; i<local_count; ++i) {
//       decomp.offsets[i] = gids_h(i)-global_min;
//     }

//     info.emplace(key,std::move(decomp));
//   }

//   return info;
// }

std::map<std::string,int>
collect_dims (const std::vector<Field>& fields,
              const std::map<decomp_key_t,DecompInfo>& decomps)
{
  std::map<std::string,int> dims;

  for (const auto& field : fields) {
    const auto& id = field.get_header().get_identifier();
    const auto& layout = id.get_layout();
    const auto& gname = id.get_grid_name();

    for (int idim=0; idim<layout.rank(); ++idim) {
      const auto name = layout.name(idim);
      const auto tag = layout.tag(idim);
      const auto dit = decomps.find(std::make_pair(gname,tag));
      const auto dimlen = dit==decomps.end() ? layout.dim(idim) : dit->second.global_extent;

      auto [it,inserted] = dims.emplace(name,dimlen);
      EKAT_REQUIRE_MSG (inserted or it->second==dimlen,
          "Error! Conflicting dimension extents while writing fields.\n"
          " - dimension name: " + name + "\n"
          " - first extent  : " + std::to_string(it->second) + "\n"
          " - new extent    : " + std::to_string(dimlen) + "\n");
    }
  }

  return dims;
}

void
define_fields (const std::string& filename,
               const std::vector<Field>& fields,
               const std::map<decomp_key_t,DecompInfo>& decomps)
{
  (void) decomps;
  DefaultMetadata meta;

  for (const auto& field : fields) {
    const auto& id = field.get_header().get_identifier();
    const auto& layout = id.get_layout();

    std::vector<std::string> dimnames;
    dimnames.reserve(layout.rank());

    for (int idim=0; idim<layout.rank(); ++idim) {
      const auto name = layout.name(idim);
      dimnames.push_back(name);
    }

    scorpio::define_var(filename,field.name(),id.get_units().to_string(),dimnames,
                        get_dtype_string(field.data_type()),
                        get_dtype_string(field.data_type()),
                        false);

    if (field.data_type()==DataType::FloatType ||
        (field.data_type()==DataType::RealType && std::is_same<Real,float>::value)) {
      scorpio::set_attribute(filename,field.name(),"_FillValue",constants::fill_value<float>);
    } else if (field.data_type()==DataType::DoubleType ||
               (field.data_type()==DataType::RealType && std::is_same<Real,double>::value)) {
      scorpio::set_attribute(filename,field.name(),"_FillValue",constants::fill_value<double>);
    }

    const auto& children = field.get_header().get_children();
    if (children.size()>0) {
      std::string children_list = "[ ";
      for (const auto& ch_w : children) {
        children_list += ch_w.lock()->get_identifier().name() + ", ";
      }
      children_list.pop_back();
      children_list.pop_back();
      children_list += " ]";
      scorpio::set_attribute(filename,field.name(),"sub_fields",children_list);
    }

    if (field.get_header().has_extra_data("io: string attributes")) {
      const auto& str_atts = field.get_header().get_extra_data<stratts_t>("io: string attributes");
      for (const auto& [att_name,att_val] : str_atts) {
        scorpio::set_attribute(filename,field.name(),att_name,att_val);
      }
      if (str_atts.count("long_name")==0) {
        scorpio::set_attribute(filename,field.name(),"long_name",meta.get_longname(field.name()));
      }
      if (str_atts.count("standard_name")==0) {
        scorpio::set_attribute(filename,field.name(),"standard_name",meta.get_standardname(field.name()));
      }
    } else {
      scorpio::set_attribute(filename,field.name(),"long_name",meta.get_longname(field.name()));
      scorpio::set_attribute(filename,field.name(),"standard_name",meta.get_standardname(field.name()));
    }
  }
}

} // anonymous namespace

namespace scream {

void
write_fields (const std::string& filename,
              const std::vector<Field>& fields)
{
  const std::vector<Field> gids_fields;
  return write_fields(filename,fields,ekat::Comm(MPI_COMM_SELF),gids_fields);
}

void
write_fields (const std::string& filename,
              const std::vector<Field>& fields,
              const ekat::Comm& comm,
              const Field& gids_field)
{
  return write_fields(filename,fields,comm,std::vector<Field>{gids_field});
}

void
write_fields (const std::string& filename,
              const FieldGroup& group)
{
  return write_fields(filename,group_to_fields(group));
}

void
write_fields (const std::string& filename,
              const FieldGroup& group,
              const ekat::Comm& comm,
              const std::vector<Field>& gids_fields)
{
  return write_fields(filename,group_to_fields(group),comm,gids_fields);
}

void
write_fields (const std::string& filename,
              const FieldGroup& group,
              const ekat::Comm& comm,
              const Field& gids_field)
{
  return write_fields(filename,group_to_fields(group),comm,gids_field);
}

void
write_fields (const std::string& filename,
              const std::vector<Field>& fields,
              const ekat::Comm& comm,
              const std::vector<Field>& decomp_dims_gids)
{
  if (fields.empty()) {
    return;
  }

  std::set<std::pair<std::string,std::string>> field_ids;
  for (const auto& field : fields) {
    const auto key = std::make_pair(field.get_header().get_identifier().get_grid_name(),field.name());
    EKAT_REQUIRE_MSG (field_ids.insert(key).second,
        "Error! Duplicate field passed to write_fields.\n"
        " - grid : " + key.first + "\n"
        " - name : " + key.second + "\n");
  }

  scorpio::register_file(filename,scorpio::Write);
  for (const auto& gids : decomp_dims_gids) {
    // Just to be safe, skip non-allocated fields
    if (not gids.is_allocated())
      continue;
    const auto& layout = gids.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG (gids.rank()==1,
        "[write_fields] Error! GIDs field must have rank 1.\n"
        " - field name: " + gids.name() + "\n"
        " - field layout: " + layout.to_string() + "\n");
    EKAT_REQUIRE_MSG (gids.data_type()==DataType::IntType,
        "[write_fields] Error! GIDs field must have data type IntType.\n"
        " - field name : " + gids.name() + "\n"
        " - field dtype: " + e2str(gids.data_type()) + "\n");

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

  const auto decomps = create_decomp_map(&comm,gids_fields);

  // Register dimensions
  for (const auto& field : fields) {
    const auto& id = field.get_header().get_identifier();
    const auto& layout = id.get_layout();
    const auto& gname = id.get_grid_name();

    for (int idim=0; idim<layout.rank(); ++idim) {
      const auto name = layout.name(idim);
      const auto tag = layout.tag(idim);
      const auto dit = decomps.find(std::make_pair(gname,tag));
      const auto dimlen = dit==decomps.end() ? layout.dim(idim) : dit->second.global_extent;

      auto [it,inserted] = dims.emplace(name,dimlen);
      EKAT_REQUIRE_MSG (inserted or it->second==dimlen,
          "Error! Conflicting dimension extents while writing fields.\n"
          " - dimension name: " + name + "\n"
          " - first extent  : " + std::to_string(it->second) + "\n"
          " - new extent    : " + std::to_string(dimlen) + "\n");
    }
  }


  for (const auto& [name,len] : dims) {
    scorpio::define_dim(filename,name,len);
  }

  define_fields(filename,fields,decomps);

  for (const auto& [key,decomp] : decomps) {
    scorpio::set_dim_decomp(filename,decomp.dim_name,decomp.offsets);
  }

  scorpio::enddef(filename);

  for (const auto& field : fields) {
    write_field_data(filename,field);
  }

  scorpio::release_file(filename);
}

} // namespace scream
