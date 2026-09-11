#include "share/data_managers/model_init.hpp"

#include <ekat_std_utils.hpp>

#include <ranges>

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params, const ekat::Comm& comm)
 : m_params (params)
 , m_comm   (comm)
{
  // Nothing to do here
}

void ModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  EKAT_REQUIRE_MSG (run_type!=RunType::Restart,
      "Error! ModelInit does not yet support restart runs.\n");

  auto gm = fm->get_grids_manager();
  for (const auto& gn : gm->get_grid_names()) {
    if (not fm->has_group("STARTUP",gn)) {
      continue;
    }

    for (auto& f : get_fields(fm,"STARTUP",gn)) {
      const auto& name = f.name();
      const auto& src_name = get_copy_source(name);
      if (src_name!="") {
        auto f_src = fm->get_field(src_name,gn);
        f.deep_copy(f_src);
        f.get_header().get_tracking().update_time_stamp(t0);
      } else if (m_params.isParameter(name)) {
        set_constant_field(f,t0);
      }
      // else: f must be inited from a startup file. Not yet supported.
    }
  }
}

std::vector<Field>
ModelInit::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name)
{
  auto group = fm->get_field_group(group_name,grid_name);
  std::vector<Field> fields;
  for (const auto& f : std::views::values(group.individual_fields())) {
    // If for some reason the field was already inited, skip it
    if (f.get_header().get_tracking().get_time_stamp().is_valid())
      continue;
    fields.push_back(f);
  }
  return fields;
}

bool ModelInit::
set_constant_field (Field& f, const util::TimeStamp& t0)
{
  const auto& name = f.name();
  if (not m_params.isParameter(name)) {
    return false;
  }

  const auto& layout = f.get_header().get_identifier().get_layout();

  // For vector fields, we allow either single value init or vector value init.
  // That is, both these are ok
  //   fname: val
  //   fname: [val1,...,valN]
  // In the first case, all entries of the field are inited to val, while in the latter,
  // each component is inited to the corresponding entry of the array.
  if (layout.is_vector_layout() and m_params.isType<std::vector<double>>(name)) {
    const auto idim = layout.get_vector_component_idx();
    const auto vec_dim = layout.get_vector_dim();
    const auto& values = m_params.get<std::vector<double>>(name);
    EKAT_REQUIRE_MSG (values.size()==static_cast<size_t>(vec_dim),
        "Error! Initial condition values array for '" + name + "' has the wrong dimension.\n"
        "       Field dimension: " + std::to_string(vec_dim) + "\n"
        "       Array dimenions: " + std::to_string(values.size()) + "\n");

    if (layout.rank()==2 && idim==1) {
      // We cannot use 'get_component' for views of rank 2 with vector dimension
      // striding fastest, since we would not get a LayoutRight view. For these views,
      // simply do a manual loop
      using kt = Field::kt_dev;
      typename kt::view_1d<double> data("data",vec_dim);
      auto data_h = Kokkos::create_mirror_view(data);
      for (int i=0; i<vec_dim; ++i) {
        data_h(i) = values[i];
      }
      Kokkos::deep_copy(data,data_h);

      const int n = layout.dim(0);
      auto v = f.get_view<double**>();
      Kokkos::parallel_for(typename kt::RangePolicy(0,n),
                           KOKKOS_LAMBDA(const int i) {
        for (int j=0; j<vec_dim; ++j) {
          v(i,j) = data(j);
        }
      });
    } else {
      // Extract a subfield for each component. This is not "too" expensive, especially
      // considering that this code is executed during initialization only.
      for (int comp=0; comp<vec_dim; ++comp) {
        auto f_i = f.get_component(comp);
        f_i.deep_copy(values[comp]);
      }
    }
  } else if (m_params.isType<int>(name)) {
    f.deep_copy(m_params.get<int>(name));
  } else {
    f.deep_copy(m_params.get<double>(name));
  }

  f.get_header().get_tracking().update_time_stamp(t0);
  return true;
}

std::string ModelInit::
get_copy_source (const std::string& name) const
{
  if (m_params.isParameter(name) and m_params.isType<std::string>(name)) {
    return m_params.get<std::string>(name);
  }
  return "";
}

} // namespace scream
