#include "share/data_managers/model_init.hpp"

#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/physics_constants.hpp"

#include <ekat_std_utils.hpp>

#include <cstdlib>
#include <ctime>
#include <functional>
#include <ranges>

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params)
 : m_params (params)
{
  // Nothing to do here
}

void ModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  auto gm = fm->get_grids_manager();
  for (const auto& gn : gm->get_grid_names()) {
    auto grid = gm->get_grid(gn);

    if (run_type==RunType::Restart) {
      if (fm->has_group("RESTART",gn)) {
        init_restart_fields(fm,grid,t0);
      }
    } else {
      // Restart files always contain the full model state, so there is no
      // need (nor a guarantee that 'topography_filename' was even given) to
      // separately load topography on a restart run.
      if (fm->has_group("TOPOGRAPHY",gn)) {
        init_topography_fields(fm,grid,t0);
      }
      if (fm->has_group("STARTUP",gn)) {
        init_startup_fields(fm,grid,t0);
      }
    }
  }

  // A restarted field is, by construction, already "perturbed" (it is
  // whatever it was at the end of the previous run), so perturbation only
  // applies to a startup run.
  if (run_type!=RunType::Restart) {
    perturb_fields(fm);
  }
}

void ModelInit::
init_startup_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::shared_ptr<const AbstractGrid>& grid,
                     const util::TimeStamp& t0)
{
  const auto& gn = grid->name();
  auto fields = get_leaf_fields(fm,"STARTUP",gn);

  // 1) Constant fields, and fields to be copied from another field (the
  //    latter is postponed, since the source field may itself need to be
  //    inited, e.g., from file, first).
  std::vector<Field> to_read;
  strmap_t<std::string> to_copy;
  for (auto& f : fields) {
    const auto& src_name = get_copy_source(f.name());
    if (src_name!="") {
      to_copy[f.name()] = src_name;
    } else if (not set_constant_field(f,t0)) {
      to_read.push_back(f);
    }
  }

  // 2) Whatever is left must come from the input file.
  if (to_read.size()>0) {
    EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
        "Error! Some fields need to be loaded from the startup file, but no "
        "'filename' was found in the input parameters.\n"
        " - grid name: " + gn + "\n");
    const auto& filename = m_params.get<std::string>("filename");
    read_fields(filename,to_read,grid,t0,get_tag_rename("STARTUP",gn));
  }

  // 3) Fields to be copied from another (by now, necessarily inited) field.
  for (const auto& [tgt_name,src_name] : to_copy) {
    auto f_tgt = fm->get_field(tgt_name,gn);
    auto f_src = fm->get_field(src_name,gn);
    f_tgt.deep_copy(f_src);
    f_tgt.get_header().get_tracking().update_time_stamp(t0);
  }

  // Some fields we inited above may be leaves of a parent field (e.g., a
  // group's monolithic field) whose time stamp does not auto-update when
  // its children's does. If all of a parent's children are now inited,
  // propagate the time stamp to the parent too.
  fixup_parents_time_stamp(fm,"STARTUP",gn,t0);
}

void ModelInit::
init_restart_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::shared_ptr<const AbstractGrid>& grid,
                     const util::TimeStamp& t0)
{
  const auto& gn = grid->name();
  auto fields = get_fields(fm,"RESTART",gn);
  if (fields.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
      "Error! Some fields need to be loaded from the restart file, but no "
      "'filename' was found in the input parameters.\n"
      " - grid name: " + gn + "\n");
  const auto& filename = m_params.get<std::string>("filename");
  read_fields(filename,fields,grid,t0,get_tag_rename("RESTART",gn));
}

void ModelInit::
init_topography_fields (const std::shared_ptr<FieldManager>& fm,
                        const std::shared_ptr<const AbstractGrid>& grid,
                        const util::TimeStamp& t0)
{
  const auto& gn = grid->name();
  auto fields = get_fields(fm,"TOPOGRAPHY",gn);
  if (fields.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("topography_filename"),
      "Error! Topography data was requested, but no 'topography_filename' "
      "was found in the input parameters.\n"
      " - grid name: " + gn + "\n");
  const auto& filename = m_params.get<std::string>("topography_filename");

  // The topography file uses different names for these fields than the
  // ones used internally by eamxx.
  const auto file_names = get_topography_file_names();
  std::vector<Field> aliased;
  for (const auto& f : fields) {
    auto it = file_names.find(f.name());
    EKAT_REQUIRE_MSG (it!=file_names.end(),
        "Error! Unrecognized topography field '" + f.name() + "'.\n");
    aliased.push_back(f.alias(it->second));
  }

  // aliased shares tracking with fields, so reading (and stamping) it also
  // stamps the fields as known by the rest of the field manager.
  read_fields(filename,aliased,grid,t0,get_tag_rename("TOPOGRAPHY",gn));
}

void ModelInit::
fixup_parents_time_stamp (const std::shared_ptr<FieldManager>& fm,
                          const std::string& group_name,
                          const std::string& grid_name,
                          const util::TimeStamp& t0)
{
  auto group = fm->get_field_group(group_name,grid_name);
  for (auto& f : std::views::values(group.individual_fields())) {
    auto& track = f.get_header().get_tracking();
    if (track.get_time_stamp().is_valid()) {
      continue;
    }
    const auto& children = track.get_children();
    if (children.size()==0) {
      continue;
    }

    bool all_inited = true;
    for (const auto& wp : children) {
      auto c = wp.lock();
      EKAT_REQUIRE_MSG (c, "Error! A weak pointer of a child field expired.\n");
      if (not c->get_time_stamp().is_valid()) {
        all_inited = false;
        break;
      }
    }
    if (all_inited) {
      track.update_time_stamp(t0);
    }
  }
}

void ModelInit::
perturb_fields (const std::shared_ptr<FieldManager>& fm)
{
  const auto perturbed_fields = m_params.get<std::vector<std::string>>("perturbed_fields",{});
  if (perturbed_fields.size()==0) {
    return;
  }

  auto gm = fm->get_grids_manager();
  EKAT_REQUIRE_MSG (gm->get_grid_names().count("physics_gll")>0,
      "Error! Random IC perturbation can only be applied to fields on the "
      "GLL grid, but no physics_gll grid was defined in the field manager.\n");
  const auto gll_grid = gm->get_grid("physics_gll");

  const auto seed = get_perturbation_seed();
  // Defines a range [1-perturbation_limit, 1+perturbation_limit] for which
  // the perturbation value will be randomly generated from.
  const auto perturbation_limit = m_params.get<Real>("perturbation_limit",0.001);
  const auto pressure_mask = build_perturbation_level_mask(gll_grid);

  const auto& gll_grid_name = gll_grid->name();
  auto dofs_gids = gll_grid->get_dofs_gids();
  for (const auto& fname : perturbed_fields) {
    auto field = fm->get_field(fname,gll_grid_name);
    EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Attempting to apply perturbation to a field that was not initialized.\n"
        "  - Field: " + fname + "\n"
        "  - Grid:  " + gll_grid_name + "\n");

    perturb(field,perturbation_limit,seed,pressure_mask,dofs_gids);
  }
}

int ModelInit::
get_perturbation_seed ()
{
  // There are two relevant params: generate_perturbation_random_seed and
  // perturbation_random_seed. We have 3 cases:
  //   1. Parameter generate_perturbation_random_seed is set true, assert perturbation_random_seed
  //      is not given and generate a random seed using std::rand() to get an integer random value.
  //   2. Parameter perturbation_random_seed is given, use this value for the seed.
  //   3. Parameter perturbation_random_seed is not given and generate_perturbation_random_seed is
  //      not given, use 0 as the random seed.
  // Case 3 is considered the default (using seed=0).
  if (m_params.get<bool>("generate_perturbation_random_seed",false)) {
    EKAT_REQUIRE_MSG (not m_params.isParameter("perturbation_random_seed"),
        "Error! Param generate_perturbation_random_seed=true, and a "
        "perturbation_random_seed is given. Only one of these can be "
        "defined for a simulation.\n");
    std::srand(std::time(nullptr));
    return std::rand();
  }
  return m_params.get<int>("perturbation_random_seed",0);
}

Field ModelInit::
build_perturbation_level_mask (const std::shared_ptr<const AbstractGrid>& gll_grid)
{
  const auto hyam_h = gll_grid->get_geometry_data("hyam").get_view<const Real*,Host>();
  const auto hybm_h = gll_grid->get_geometry_data("hybm").get_view<const Real*,Host>();
  constexpr auto ps0 = physics::Constants<Real>::P0.value;
  const auto min_pressure = m_params.get<Real>("perturbation_minimum_pressure",1050.0);

  using namespace ShortFieldTagsNames;
  const auto& pmask_lt = gll_grid->get_vertical_layout(LEV);
  const auto nondim = ekat::units::none;
  FieldIdentifier pmask_fid("lev_mask",pmask_lt,nondim,gll_grid->name(),DataType::IntType);
  Field pressure_mask(pmask_fid,true);
  auto pmask_h = pressure_mask.get_view<int*,Host>();
  for (int ilev=0; ilev<pmask_lt.dim(0); ++ilev) {
    const auto pref = (hyam_h(ilev)*ps0 + hybm_h(ilev)*ps0)/100; // Reference pressure ps0 is in Pa, convert to millibar
    pmask_h(ilev) = static_cast<int>(pref > min_pressure);
  }
  pressure_mask.sync_to_dev();
  return pressure_mask;
}

std::vector<Field>
ModelInit::
get_leaf_fields (const std::shared_ptr<FieldManager>& fm,
                 const std::string& group_name,
                 const std::string& grid_name)
{
  auto group = fm->get_field_group(group_name,grid_name);
  std::vector<Field> fields;

  std::function<void(const Field&)> collect_leaves = [&] (const Field& f) {
    if (f.get_header().get_tracking().get_time_stamp().is_valid()) {
      return;
    }
    const auto& children = f.get_header().get_children();
    if (children.size()>0) {
      for (const auto& wp : children) {
        auto c = wp.lock();
        EKAT_REQUIRE_MSG (c, "Error! A weak pointer of a child field expired.\n");
        collect_leaves(fm->get_field(c->get_identifier().name(),grid_name));
      }
    } else {
      fields.push_back(f);
    }
  };
  for (const auto& f : std::views::values(group.individual_fields())) {
    collect_leaves(f);
  }
  return fields;
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
    if (f.get_header().get_tracking().get_time_stamp().is_valid()) {
      continue;
    }
    // If the parent of f is also part of this group, skip f: reading (and
    // stamping) the parent will automatically take care of f too.
    auto p = f.get_header().get_parent();
    if (p and ekat::contains(p->get_tracking().get_groups_names(),group_name)) {
      continue;
    }
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

ModelInit::strmap_t<std::string> ModelInit::
get_tag_rename (const std::string& group_name,
                const std::string& grid_name) const
{
  if (group_name=="TOPOGRAPHY" and grid_name=="physics_gll") {
    return { {"ncol","ncol_d"} };
  }
  return {};
}

ModelInit::strmap_t<std::string> ModelInit::
get_topography_file_names () const
{
  return {
    {"phis",  "PHIS_d"},
    {"sgh30", "SGH30"},
    {"sgh",   "SGH"}
  };
}

std::string ModelInit::
get_copy_source (const std::string& name) const
{
  if (m_params.isParameter(name) and m_params.isType<std::string>(name)) {
    return m_params.get<std::string>(name);
  }
  return "";
}

void ModelInit::
read_fields (const std::string& filename,
            std::vector<Field>& fields,
            const std::shared_ptr<const AbstractGrid>& grid,
            const util::TimeStamp& t0,
            const strmap_t<std::string>& tag_rename)
{
  FieldReader reader;
  reader.set_fields(fields);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),grid->get_comm());
  reader.set_file_specs(filename,tag_rename);
  reader.read();

  for (auto& f : fields) {
    f.get_header().get_tracking().update_time_stamp(t0);
  }
}

} // namespace scream
