#include "dynamics/homme/homme_model_init.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_instance.hpp"
#include "dynamics/homme/interface/eamxx_homme_interface.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"

#include <ekat_assert.hpp>
#include <ekat_std_utils.hpp>

namespace scream
{

namespace {

// Parse a name of the form "physics_pgn". Return -1 if not an FV physics grid
// name, otherwise N in pgN.
int get_phys_grid_fv_param (const std::string& grid_name)
{
  const std::string prefix = "physics_pg";
  if (grid_name.size()<=prefix.size() or grid_name.substr(0,prefix.size())!=prefix) {
    return -1;
  }
  try {
    return std::stoi(grid_name.substr(prefix.size()));
  } catch (...) {
    return -1;
  }
}

} // anonymous namespace

void HommeModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  auto instance = HommeInstance::get();
  if (not instance->has_dyn_field("v_dyn")) {
    // The Homme atm process is not in use (e.g., only the grids are), so there is no Homme state to init.
    ModelInit::run(fm,t0,run_type);
    return;
  }

  auto gm = fm->get_grids_manager();
  const std::string gll = "physics_gll";

  // Is the physics grid a FV one?
  m_pgN = -1;
  m_pg_grid = "";
  for (const auto& gn : gm->get_grid_names()) {
    const int N = get_phys_grid_fv_param(gn);
    if (N>0) {
      m_pgN = N;
      m_pg_grid = gn;
    }
  }
  const bool pg = m_pgN>0;
  const bool restart = run_type==RunType::Restart;

  const auto& params = Homme::Context::singleton().get<Homme::SimulationParams>();

  // 1) Init the atm fields. The field manager holding the GLL ICs is the main one,
  //    unless we are running on a FV physics grid.
  auto gll_ic_fm = fm;
  if (restart) {
    init_all_grids(fm,t0,run_type);
    load_restart_extra_data();
  } else if (pg) {
    gll_ic_fm = create_gll_ic_field_manager(fm);
    init_all_grids(gll_ic_fm,t0,run_type);
    finalize_startup(gll_ic_fm,t0);

    // Init the remaining fields (e.g., constants, topography, trace gases),
    // but don't do IOP/perturbation again.
    init_all_grids(fm,t0,run_type);
  } else {
    init_all_grids(fm,t0,run_type);
    finalize_startup(fm,t0);
  }

  // 2) Complete Homme prim_init1_xyz sequence
  prim_complete_init1_phase_f90 ();

  // 3) Init the gll<->pg remapper, and set the views of the dyn fields in Homme's data structures
  if (pg) {
    instance->init_fv_phys(m_pgN);
  }
  instance->init_views ();

  // The first D->P remap during initialization happens before the dycore has computed
  // the strain components, so start from a benign value. Homme overwrites them after each
  // dynamics step when 3D turbulence is enabled.
  if (pg and params.do_3d_turbulence) {
    instance->get_dyn_field("shear_strain3d_components_dyn").deep_copy(0);
    fm->get_field("tke_shear_strain3d_components",m_pg_grid).deep_copy(0);
  }

  // 4) Init the dyn state
  if (restart) {
    // For BFB restarts, set nstep counter in Homme's TimeLevel to match the restarted value.
    EKAT_REQUIRE_MSG (m_restart_extra_data.count("homme_nsteps")==1,
        "Error! The restart extra data 'homme_nsteps' is needed to restart Homme.\n");
    const auto& nstep = std::any_cast<const int&>(*m_restart_extra_data.at("homme_nsteps"));
    Homme::Context::singleton().get<Homme::TimeLevel>().nstep = nstep;
    set_homme_param("num_steps",nstep);

    instance->finish_dyn_state_restart ();
  } else {
    // Init Homme's state from the ICs on the GLL grid
    auto gll2dyn = gm->create_remapper(gll,"dynamics");
    auto tracers_gll = gll_ic_fm->get_field_group("tracers",gll).monolithic_field();

    if (pg) {
      // The monolithic tracers field in the IC field manager may not have the tracers in the
      // same order as the one on the pg grid (and Homme's tracers follow the latter),
      // so copy the tracers (by name) in a field with the right order.
      const auto pg_tracers = fm->get_field_group("tracers",m_pg_grid);
      const int nq = pg_tracers.size();
      FieldIdentifier fid("tracers",gm->get_grid(gll)->get_3d_vector_layout(ShortFieldTagsNames::LEV,nq),
                          ekat::units::none,gll);
      Field ordered(fid);
      ordered.get_header().get_alloc_properties().request_allocation(HOMMEXX_PACK_SIZE);
      ordered.allocate_view();
      for (const auto& [name,f] : pg_tracers.individual_fields()) {
        const int idx = f.get_header().get_alloc_properties().get_subview_info().slice_idx;
        ordered.get_component(idx).deep_copy(gll_ic_fm->get_field(name,gll));
      }
      tracers_gll = ordered;
    }

    auto pseudo_density_gll = gll_ic_fm->get_field("pseudo_density",gll);
    instance->init_dyn_state_from_gll_ic(
        *gll2dyn,
        gll_ic_fm->get_field("horiz_winds",gll),
        gll_ic_fm->get_field("T_mid",gll),
        gll_ic_fm->get_field("ps",gll),
        gll_ic_fm->get_field("phis",gll),
        tracers_gll,
        pseudo_density_gll);

    // The pseudo density is computed by Homme (it is not in the IC file)
    pseudo_density_gll.get_header().get_tracking().update_time_stamp(t0);
  }

  // 5) Complete homme model initialization
  prim_init_model_f90 ();

  // 6) For FV physics grids, remap the dyn state to the physics grid
  if (pg) {
    const auto& pgn = m_pg_grid;
    auto get = [&](const std::string& name) { return fm->get_field(name,pgn); };

    Field strain3d_fv;
    if (params.do_3d_turbulence) {
      strain3d_fv = get("tke_shear_strain3d_components");
    }
    auto tracers = fm->get_field_group("tracers",pgn).monolithic_field();

    auto T_prev  = instance->get_phys_field("FT_phys");
    auto uv_prev = instance->get_phys_field("FM_phys");

    if (restart) {
      // The state fields (T_mid, horiz_winds, and tracers) of the pg grid were restarted, and
      // they must not be overwritten, since they contain data used by HommeDynamics to compute
      // tendencies. So remap into temporary copies. We still need the remapped T and uv, to
      // set the "previous step" ones.
      auto make_tmp = [](const Field& f) {
        Field tmp(f.get_header().get_identifier());
        tmp.get_header().get_alloc_properties().request_allocation(HOMMEXX_PACK_SIZE);
        tmp.allocate_view();
        return tmp;
      };
      auto T_tmp  = make_tmp(get("T_mid"));
      auto uv_tmp = make_tmp(get("horiz_winds"));
      auto q_tmp  = make_tmp(tracers);

      instance->remap_dyn_to_fv_phys(get("ps"),get("phis"),T_tmp,get("omega"),
                                     uv_tmp,q_tmp,get("pseudo_density"),strain3d_fv);
      T_prev.deep_copy(T_tmp);
      uv_prev.deep_copy(uv_tmp);
    } else {
      instance->remap_dyn_to_fv_phys(get("ps"),get("phis"),get("T_mid"),get("omega"),
                                     get("horiz_winds"),tracers,get("pseudo_density"),strain3d_fv);
      T_prev.deep_copy(get("T_mid"));
      uv_prev.deep_copy(get("horiz_winds"));

      // The remapped fields are not inited by the IC file, so we need to set their time stamp
      std::vector<Field> remapped = {
        get("T_mid"), get("horiz_winds"), get("ps"), get("phis"), get("omega"),
        get("pseudo_density"), tracers
      };
      if (params.do_3d_turbulence) {
        remapped.push_back(strain3d_fv);
      }
      update_time_stamps(remapped,t0);
    }
  }
}

std::vector<Field> HommeModelInit::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name)
{
  auto fields = ModelInit::get_fields(fm,group_name,grid_name);
  if (m_pgN<0) {
    return fields;
  }

  if (grid_name=="physics_gll" and group_name=="RESTART") {
    // physics_gll only has the trace gases ICs (see eamxx_homme_fv_phys.cpp), which are
    // not part of a restart file (the pg grid ones are).
    return {};
  }

  if (grid_name==m_pg_grid and group_name=="STARTUP") {
    // The IC file only has data on the GLL grid. The state fields are inited by remapping
    // the dyn state, so skip them. In fact, skip anything that is not constant or topography.
    const auto& overrides = m_params.sublist("constant_fields").sublist("overrides");
    const auto topo_names = get_topography_file_names(grid_name);
    std::erase_if(fields,[&](const Field& f) {
      return not (overrides.isParameter(f.name()) or topo_names.count(f.name())>0);
    });
  }
  return fields;
}

HommeModelInit::strmap_t<std::string> HommeModelInit::
get_topography_file_names (const std::string& grid_name) const
{
  auto names = ModelInit::get_topography_file_names(grid_name);
  if (grid_name==m_pg_grid) {
    names.erase("phis");
  }
  return names;
}

std::shared_ptr<FieldManager> HommeModelInit::
create_gll_ic_field_manager (const std::shared_ptr<FieldManager>& fm) const
{
  using namespace ShortFieldTagsNames;
  constexpr int N = HOMMEXX_PACK_SIZE;

  auto gm = fm->get_grids_manager();
  auto gll = gm->get_grid("physics_gll");
  const auto& gln = gll->name();
  const auto& pgn = m_pg_grid;

  auto units = [&](const std::string& name) {
    return fm->get_field(name,pgn).get_header().get_identifier().get_units();
  };

  auto ic_fm = std::make_shared<FieldManager>(gm);

  // The state fields that Homme owns on the pg grid
  ic_fm->register_field(FieldRequest(FieldIdentifier("horiz_winds",gll->get_3d_vector_layout(LEV,2),units("horiz_winds"),gln),N));
  ic_fm->register_field(FieldRequest(FieldIdentifier("T_mid",gll->get_3d_scalar_layout(LEV),units("T_mid"),gln),N));
  ic_fm->register_field(FieldRequest(FieldIdentifier("ps",gll->get_2d_scalar_layout(),units("ps"),gln)));
  ic_fm->register_field(FieldRequest(FieldIdentifier("phis",gll->get_2d_scalar_layout(),units("phis"),gln)));

  // The pseudo density is not an IC: it is computed (hydrostatically) from ps
  ic_fm->register_field(FieldRequest(FieldIdentifier("pseudo_density",gll->get_3d_scalar_layout(LEV),units("pseudo_density"),gln),N));

  // Tracers. All of them are inited from the IC file/constants.
  const auto pg_tracers = fm->get_field_group("tracers",pgn);
  std::vector<std::string> tracers_names;
  for (const auto& [name,f] : pg_tracers.individual_fields()) {
    const auto& u = f.get_header().get_identifier().get_units();
    ic_fm->register_field(FieldRequest(FieldIdentifier(name,gll->get_3d_scalar_layout(LEV),u,gln),{"tracers"},N));
    tracers_names.push_back(name);
  }
  ic_fm->register_group(GroupRequest("tracers",gln,N,MonolithicAlloc::Required));
  ic_fm->register_group(GroupRequest("STARTUP",gln));
  ic_fm->registration_ends();

  // As the driver does, add U/V subfields of horiz_winds
  auto hw = ic_fm->get_field("horiz_winds",gln);
  ic_fm->add_field(hw.get_component(0).alias("U"));
  ic_fm->add_field(hw.get_component(1).alias("V"));

  // Only leaf fields are in the STARTUP group
  for (const auto& name : {"U","V","T_mid","ps","phis"}) {
    ic_fm->add_to_group(name,gln,"STARTUP");
  }
  for (const auto& name : tracers_names) {
    ic_fm->add_to_group(name,gln,"STARTUP");
  }

  return ic_fm;
}

} // namespace scream
