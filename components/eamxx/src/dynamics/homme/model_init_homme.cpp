#include "dynamics/homme/model_init_homme.hpp"
#include "dynamics/homme/eamxx_homme_context.hpp"
#include "dynamics/homme/eamxx_homme_process_interface.hpp"

namespace scream
{

void ModelInitHomme::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  // Do the "normal" field init (STARTUP/RESTART/TOPOGRAPHY, on all grids),
  // exactly like any other atm process. The base/derived get_leaf_fields and
  // get_fields hooks take care of any Homme-specific field selection (e.g.,
  // ModelInitPG2 skips the physics_pg2-grid fields that Homme's own remap,
  // rather than the initial condition file, is responsible for).
  ModelInit::run(fm,t0,run_type);

  if (run_type==RunType::Restart) {
    // Homme's dynamics-grid state is part of the restart file (read via the
    // base ModelInit::run() above, through the "dynamics" grid's RESTART
    // group), but only its n0 time level: since Homme::Context's own views
    // are, by this point, already aliased to this data (see
    // create_requests/set_computed_group_impl), make it fully consistent
    // right away -- Homme's own initialization (prim_complete_init1_phase_f90,
    // called from HommeDynamics::initialize_impl before it gets to call
    // restart_homme_state, its own, later, otherwise-redundant version of
    // this) needs that, not just n0.
    finish_homme_dyn_state_restart(fm->get_grids_manager(),
                                    EamxxHommeContext::singleton().dyn_state_fields);
    return;
  }

  // Populate Homme's dynamics-grid state from the (now-inited) physics_gll
  // initial condition. By this point (well after HommeDynamics's
  // create_requests()/set_computed_group_impl()), the dyn_state_fields
  // EamxxHommeContext holds are already aliased to Homme::Context's own
  // ElementsState/Tracers views, so this also directly populates what
  // Homme's dynamical core will use.
  const std::string rgn = "physics_gll";
  auto pseudo_density_gll = fm->get_field("pseudo_density",rgn);
  init_homme_dyn_state_from_gll_ic(
      fm->get_grids_manager(),
      fm->get_field("horiz_winds",rgn),
      fm->get_field("T_mid",rgn),
      fm->get_field("ps",rgn),
      fm->get_field("phis",rgn),
      fm->get_field_group("tracers",rgn).monolithic_field(),
      pseudo_density_gll,
      EamxxHommeContext::singleton().dyn_state_fields,
      t0);
}

} // namespace scream
