#include "dynamics/homme/homme_model_init.hpp"
#include "dynamics/homme/homme_instance.hpp"
#include "dynamics/homme/interface/eamxx_homme_interface.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "TimeLevel.hpp"

#include <ekat_assert.hpp>

namespace scream
{

void HommeModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  // Init all the atm fields (and the dyn-grid fields in the restart group, for restart runs)
  ModelInit::run(fm,t0,run_type);

  auto instance = HommeInstance::get();
  if (not instance->has_dyn_field("v_dyn")) {
    // The Homme atm process is not in use (e.g., only the grids are), so there is no Homme state to init.
    return;
  }

  auto gm = fm->get_grids_manager();
  const std::string gll = "physics_gll";

  // Complete Homme prim_init1_xyz sequence
  prim_complete_init1_phase_f90 ();

  // Sets the scream views into the hommexx internal data structures
  instance->init_views ();

  if (run_type==RunType::Initial) {
    // Init Homme's state from the ICs on the GLL grid
    auto gll2dyn = gm->create_remapper(gll,"dynamics");
    auto pseudo_density = fm->get_field("pseudo_density",gll);
    instance->init_dyn_state_from_gll_ic(
        *gll2dyn,
        fm->get_field("horiz_winds",gll),
        fm->get_field("T_mid",gll),
        fm->get_field("ps",gll),
        fm->get_field("phis",gll),
        fm->get_field_group("tracers",gll).monolithic_field(),
        pseudo_density);

    // The pseudo density is computed by Homme (it is not in the IC file)
    pseudo_density.get_header().get_tracking().update_time_stamp(t0);
  } else {
    // For BFB restarts, set nstep counter in Homme's TimeLevel to match the restarted value.
    EKAT_REQUIRE_MSG (m_restart_extra_data.count("homme_nsteps")==1,
        "Error! The restart extra data 'homme_nsteps' is needed to restart Homme.\n");
    const auto& nstep = std::any_cast<const int&>(*m_restart_extra_data.at("homme_nsteps"));
    Homme::Context::singleton().get<Homme::TimeLevel>().nstep = nstep;
    set_homme_param("num_steps",nstep);

    instance->finish_dyn_state_restart ();
  }

  // Complete homme model initialization
  prim_init_model_f90 ();
}

} // namespace scream
