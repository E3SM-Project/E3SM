#include "share/data_managers/model_init.hpp"

#include "share/data_managers/IOPDataManager.hpp"
#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_assert.hpp>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <ranges>

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params)
 : m_params (params)
{
  // The constants are in the 'constant_fields' sublist, which has two sublists
  // (if not present, they are created empty): 'overrides' and 'defaults'.
  m_params.sublist("constant_fields").sublist("overrides");
  m_params.sublist("constant_fields").sublist("defaults");
}

void ModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  init_all_grids(fm,t0,run_type);

  if (run_type==RunType::Restart) {
    load_restart_extra_data();
  } else {
    finalize_startup(fm,t0);
  }
}

void ModelInit::
init_all_grids (const std::shared_ptr<FieldManager>& fm,
                const util::TimeStamp& t0,
                const RunType run_type)
{
  const bool restart = run_type==RunType::Restart;

  if (not restart) {
    check_no_constants_for_composite_fields(fm);
  }

  const auto group_name = restart ? "RESTART" : "STARTUP";
  for (const auto& gn : fm->get_grids_manager()->get_grid_names()) {
    if (not fm->has_group(group_name,gn)) {
      continue;
    }
    if (restart) {
      init_restart_fields(fm,gn,t0);
    } else {
      init_startup_fields(fm,gn,t0);
    }
  }
}

void ModelInit::
finalize_startup (const std::shared_ptr<FieldManager>& fm,
                  const util::TimeStamp& t0)
{
  if (m_iop_data_manager) {
    // Load IOP data file data for initial time stamp
    m_iop_data_manager->read_iop_file_data(t0);

    // Now that ICs are processed, set appropriate fields using IOP file data.
    // Since ICs are loaded on the GLL grid, we set those fields only, and dynamics
    // will take care of the rest (for PG2 case).
    if (fm->get_grids_manager()->get_grid_names().count("physics_gll")>0) {
      m_iop_data_manager->set_fields_from_iop_data(fm,"physics_gll");
    }
  }

  perturb_fields(fm);
}

void ModelInit::
init_startup_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::string& grid_name,
                     const util::TimeStamp& t0)
{
  const auto& constants_pl = m_params.sublist("constant_fields");
  const auto& overrides_pl = constants_pl.sublist("overrides");
  const auto& defaults_pl  = constants_pl.sublist("defaults");

  // Check if a field is in the IC file (if any). The file is opened on the fly (and released
  // right away), and only if we need to check for a field with a default constant.
  const bool has_ic_file = m_params.isParameter("filename");
  const std::string ic_file = has_ic_file ? m_params.get<std::string>("filename") : "";
  auto in_ic_file = [&](const std::string& fname) {
    return has_ic_file and scorpio::has_var(ic_file,fname);
  };

  // Process all fields in the STARTUP group. For each, either init to
  // a constant (if provided), or add it to the list of fields to read from
  // the IC file (or from the topography file, for topography fields).
  // NOTE: default constants are not used for the topography fields.
  const auto topo_names = get_topography_file_names(grid_name);
  std::vector<Field> inited;
  std::vector<Field> to_read;
  std::vector<Field> topo_to_read;
  for (auto f : get_fields(fm,"STARTUP",grid_name)) {
    const auto& fname = f.name();
    if (overrides_pl.isParameter(fname)) {
      set_constant_field(f,overrides_pl);
      inited.push_back(f);
    } else if (topo_names.count(fname)>0) {
      // The field has a different name in the topography file. The alias shares
      // the data with f, so reading into the alias inits f.
      topo_to_read.push_back(f.alias(topo_names.at(fname)));
      inited.push_back(f);
    } else if (defaults_pl.isParameter(fname) and not in_ic_file(fname)) {
      set_constant_field(f,defaults_pl);
      inited.push_back(f);
    } else {
      to_read.push_back(f);
    }
  }

  if (to_read.size()>0) {
    EKAT_REQUIRE_MSG (has_ic_file,
        "Error! Some fields need to be loaded from the IC file, but no 'filename' "
        "was found in the initial conditions parameters.\n"
        " - grid name: " + grid_name + "\n"
        " - fields: " + to_read.front().name() + (to_read.size()>1 ? ", ..." : "") + "\n");
    read_fields(ic_file,to_read,fm->get_grids_manager()->get_grid(grid_name));
    inited.insert(inited.end(),to_read.begin(),to_read.end());
  }

  if (topo_to_read.size()>0) {
    EKAT_REQUIRE_MSG (m_params.isParameter("topography_filename"),
        "Error! Topography data was requested in the FieldManager, but no "
        "'topography_filename' was found in the initial conditions parameters "
        "(nor was a constant provided for the field).\n"
        " - grid name: " + grid_name + "\n"
        " - field: " + topo_to_read.front().name() + "\n");
    read_fields(m_params.get<std::string>("topography_filename"),topo_to_read,
                fm->get_grids_manager()->get_grid(grid_name),
                get_topography_tag_rename(grid_name));
  }

  update_time_stamps(inited,t0);
}

void ModelInit::
init_restart_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::string& grid_name,
                     const util::TimeStamp& t0)
{
  auto fields = get_fields(fm,"RESTART",grid_name);
  if (fields.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
      "Error! Some fields need to be loaded from the restart file, but no 'filename' "
      "was found in the initial conditions parameters.\n"
      " - grid name: " + grid_name + "\n");
  const auto& filename = m_params.get<std::string>("filename");

  try {
    // IOP only affects startup runs: a restart file has the full state of every column
    read_fields(filename,fields,fm->get_grids_manager()->get_grid(grid_name),{},false);
  } catch (std::exception& e) {
    EKAT_ERROR_MSG (
        std::string(e.what()) + "\n"
        "NOTE: if the error above complains about a missing variable 'U' or 'V',\n"
        "      you may be trying to restart from a restart file created before\n"
        "      PR #8606 (https://github.com/E3SM-Project/E3SM/pull/8606), which\n"
        "      changed EAMxx to init/restart 'U'/'V' directly, rather than the\n"
        "      vector field 'horiz_winds'. See the 'Restart file compatibility'\n"
        "      section in components/eamxx/docs/user/model_configuration.md for\n"
        "      instructions on how to convert an old restart file.\n");
  }

  update_time_stamps(fields,t0);
}

void ModelInit::
load_restart_extra_data () const
{
  if (m_restart_extra_data.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
      "Error! Restart extra data needs to be loaded from the restart file, but no 'filename' "
      "was found in the initial conditions parameters.\n");
  const auto& filename = m_params.get<std::string>("filename");

  for (const auto& it : m_restart_extra_data) {
    const auto& name = it.first;
          auto& any  = *it.second;

    if (any.type()==typeid(int)) {
      std::any_cast<int&>(any) = scorpio::get_attribute<int>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(std::int64_t)) {
      std::any_cast<std::int64_t&>(any) = scorpio::get_attribute<std::int64_t>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(float)) {
      std::any_cast<float&>(any) = scorpio::get_attribute<float>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(double)) {
      std::any_cast<double&>(any) = scorpio::get_attribute<double>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(std::string)) {
      std::any_cast<std::string&>(any) = scorpio::get_attribute<std::string>(filename,"GLOBAL",name);
    } else {
      EKAT_ERROR_MSG (
          "Error! Unrecognized/unsupported concrete type for restart extra data.\n"
          " - extra data name  : " + name + "\n"
          " - extra data typeid: " + std::string(any.type().name()) + "\n");
    }
  }
}

std::vector<Field> ModelInit::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name)
{
  std::vector<Field> fields;
  auto group = fm->get_field_group(group_name,grid_name);
  for (const auto& f : std::views::values(group.individual_fields())) {
    // Skip fields that were already inited
    if (not f.get_header().get_tracking().get_time_stamp().is_valid()) {
      fields.push_back(f);
    }
  }
  return fields;
}

void ModelInit::
read_fields (const std::string& filename,
             const std::vector<Field>& fields,
             const std::shared_ptr<const AbstractGrid>& grid,
             const strmap_t<std::string>& tag_rename,
             const bool use_iop)
{
  if (m_iop_data_manager and use_iop) {
    // For IOP enabled, we load from file and copy data from the closest
    // lat/lon column to every other column.
    // NOTE: the topography file only contains lat/lon on the PG2 grid, while
    //       the IC file contains them on the GLL grid. So, if we need to
    //       read a topography field on the GLL grid, we set up the io info
    //       using the IC file.
    const auto& gn = grid->name();
    if (m_iop_io_grids.count(gn)==0) {
      const auto& info_file = m_params.get<std::string>(gn=="physics_gll" ? "filename" : "topography_filename");
      m_iop_data_manager->setup_io_info(info_file,grid);
      m_iop_io_grids.insert(gn);
    }
    m_iop_data_manager->read_fields_from_file_for_iop(filename,fields,grid);
    return;
  }

  FieldReader reader;
  reader.set_file_specs(filename,tag_rename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),grid->get_comm());
  reader.set_fields(fields);
  reader.read();
}

void ModelInit::
update_time_stamps (const std::vector<Field>& fields,
                    const util::TimeStamp& t0) const
{
  for (auto f : fields) {
    f.get_header().get_tracking().update_time_stamp(t0);

    // Also init the timestamp of children (if any)
    for (auto& c : f.get_header().get_children()) {
      c.lock()->get_tracking().update_time_stamp(t0);
    }

    // If there is a parent, and all children were inited,
    // also init the parent's timestamp
    if (auto p = f.get_header().get_parent(); p!=nullptr) {
      bool p_inited = true;
      for (auto c : p->get_children()) {
        p_inited &= c.lock()->get_tracking().get_time_stamp().is_valid();
      }
      if (p_inited) {
        p->get_tracking().update_time_stamp(t0);
      }
    }
  }
}

ModelInit::strmap_t<std::string> ModelInit::
get_topography_file_names (const std::string& /* grid_name */) const
{
  return {
    {"phis",  "PHIS_d"},
    {"sgh30", "SGH30"},
    {"sgh",   "SGH"}
  };
}

ModelInit::strmap_t<std::string> ModelInit::
get_topography_tag_rename (const std::string& grid_name) const
{
  if (grid_name=="physics_gll") {
    return { {"ncol","ncol_d"} };
  }
  return {};
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
      "Error! Random perturbation can only be applied to fields on the GLL grid, "
      "but no physics_gll grid was defined in the FieldManager.\n");
  const auto gll_grid = gm->get_grid("physics_gll");
  const auto& gll_grid_name = gll_grid->name();

  const auto seed = get_perturbation_seed();

  // Defines a range [1-perturbation_limit, 1+perturbation_limit] for which
  // the perturbation value will be randomly generated from.
  const auto perturbation_limit = m_params.get<Real>("perturbation_limit",0.001);
  const auto pressure_mask = build_perturbation_level_mask(gll_grid);

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
get_perturbation_seed () const
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
        "Error! Param generate_perturbation_random_seed=true, and "
        "a perturbation_random_seed is given. Only one of these can "
        "be defined for a simulation.\n");
    std::srand(std::time(nullptr));
    return std::rand();
  }
  return m_params.get<int>("perturbation_random_seed",0);
}

Field ModelInit::
build_perturbation_level_mask (const std::shared_ptr<const AbstractGrid>& gll_grid) const
{
  // Define a level mask using reference pressure and the perturbation_minimum_pressure parameter.
  // This mask dictates which levels we apply a perturbation.
  const auto hyam_h = gll_grid->get_geometry_data("hyam").get_view<const Real*,Host>();
  const auto hybm_h = gll_grid->get_geometry_data("hybm").get_view<const Real*,Host>();
  constexpr auto ps0 = physics::Constants<Real>::P0.value;
  const auto min_pressure = m_params.get<Real>("perturbation_minimum_pressure",1050.0);

  using namespace ShortFieldTagsNames;
  const auto& pmask_lt = gll_grid->get_vertical_layout(LEV);
  FieldIdentifier pmask_fid("lev_mask",pmask_lt,ekat::units::none,gll_grid->name(),DataType::IntType);
  Field pressure_mask(pmask_fid,true);
  auto pmask_h = pressure_mask.get_view<int*,Host>();
  for (int ilev=0; ilev<pmask_lt.dim(0); ++ilev) {
    // Reference pressure ps0 is in Pa, convert to millibar
    const auto pref = (hyam_h(ilev)*ps0 + hybm_h(ilev)*ps0)/100;
    pmask_h(ilev) = static_cast<int>(pref > min_pressure);
  }
  pressure_mask.sync_to_dev();
  return pressure_mask;
}

void ModelInit::
check_no_constants_for_composite_fields (const std::shared_ptr<FieldManager>& fm)
{
  const auto& constants_pl = m_params.sublist("constant_fields");
  const auto& overrides_pl = constants_pl.sublist("overrides");
  const auto& defaults_pl  = constants_pl.sublist("defaults");

  for (const auto& gn : fm->get_grids_manager()->get_grid_names()) {
    for (const auto& it : fm->get_repo(gn)) {
      const auto& f = *it.second;
      const auto& children = f.get_header().get_children();
      if (children.size()>0 and
          (overrides_pl.isParameter(f.name()) or defaults_pl.isParameter(f.name()))) {
        std::string child_names;
        for (auto c : children) {
          child_names += c.lock()->get_identifier().name() + " ";
        }
        EKAT_ERROR_MSG (
            "Error! Cannot set an initial condition for field '" + f.name() + "' directly, "
            "since it has subfields.\n"
            "       Grid name:     " + gn + "\n"
            "       Field:         " + f.name() + "\n"
            "       Subfields:     " + child_names + "\n"
            "       Please, set the initial condition for each subfield individually.\n");
      }
    }
  }
}

void ModelInit::
set_constant_field (Field& f, const ekat::ParameterList& constants_pl) const
{
  const auto& fid = f.get_header().get_identifier();
  const auto& layout = fid.get_layout();
  const auto& name = fid.name();

  // The constant can be specified in one of these forms
  //   fname: val
  //   fname: [val]
  //   fname: [val1,...,valN]
  // where the values are reals (write 1.0, not 1, since the ekat yaml parser reads 1 as an int).
  // The XML defaults enforce that all constants are array(real), so that is what atmchange
  // and buildnml produce. The scalar form is only accepted to keep the input.yaml files of
  // standalone tests simple(r).
  // In the first two cases, all entries of the field are inited to val (for any field).
  // In the last case, the field must be a vector field with N components, and each
  // component is inited to the corresponding entry.
  std::vector<double> values;
  if (constants_pl.isType<std::vector<double>>(name)) {
    values = constants_pl.get<std::vector<double>>(name);
  } else if (constants_pl.isType<double>(name)) {
    values.push_back(constants_pl.get<double>(name));
  } else {
    EKAT_ERROR_MSG ("Error! Invalid initial condition for field '" + name + "': only scalar "
                    "or vector double arguments are allowed.\n"
                    "Note: If you used int values, convert to double (e.g., write 1.0 instead of 1).\n");
  }

  // Check all inputs before modifying the field
  const int vec_dim = layout.is_vector_layout() ? layout.get_vector_dim() : 1;
  EKAT_REQUIRE_MSG (values.size()==1 or values.size()==static_cast<size_t>(vec_dim),
      "Error! Initial condition values array for '" + name + "' has the wrong dimension.\n"
      "       Grid name:       " + fid.get_grid_name() + "\n"
      "       Field layout:    " + layout.to_string() + "\n"
      "       Array dimension: " + std::to_string(values.size()) + "\n"
      "       Valid dimensions: 1 (set all entries to the same value)" +
      (vec_dim>1 ? " or " + std::to_string(vec_dim) + " (one value per vector component)\n" : "\n"));

  // Integer fields can be set from (real) values, but the conversion must be exact
  // (e.g., 2.0 is ok, 2.5 throws), since we do not want to silently truncate. The value
  // must also fit in an int (this also rules out inf, and NaN fails the first check).
  const bool is_int = f.data_type()==DataType::IntType;
  if (is_int) {
    for (const auto v : values) {
      EKAT_REQUIRE_MSG (v==std::floor(v) and std::abs(v)<=std::numeric_limits<int>::max(),
          "Error! Cannot initialize the integer field '" + name + "' with a non-integral value, "
          "or a value outside the range of int.\n"
          "       Value: " + std::to_string(v) + "\n");
    }
  }

  auto set_value = [&](const Field& x, const double v) {
    if (is_int) {
      x.deep_copy(static_cast<int>(v));
    } else {
      x.deep_copy(v);
    }
  };

  if (values.size()==1) {
    set_value(f,values[0]);
  } else {
    // Extract a subfield for each component. This is not "too" expensive, especially
    // considering that this code is executed during initialization only.
    for (int comp=0; comp<vec_dim; ++comp) {
      set_value(f.get_component(comp),values[comp]);
    }
  }
}

} // namespace scream
