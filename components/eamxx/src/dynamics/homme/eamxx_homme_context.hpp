#ifndef EAMXX_HOMME_CONTEXT_HPP
#define EAMXX_HOMME_CONTEXT_HPP

#include "share/field/field.hpp"

#include <map>
#include <string>

namespace scream
{

// A small, Homme-package-scoped registry (in spirit similar to Homme::Context,
// but for state that bridges eamxx and Homme) used to share HommeDynamics's
// dynamics-grid state Field objects with ModelInitHomme (and its subclasses),
// without coupling ModelInit to the HommeDynamics atmosphere process itself.
//
// HommeDynamics populates dyn_state_fields, during create_requests(), with
// the same (shared-header) Field objects it uses internally for the
// dynamics-grid prognostic state ("v_dyn", "vtheta_dp_dyn", "dp3d_dyn",
// "w_int_dyn", "phi_int_dyn", "ps_dyn", "phis_dyn", "omega_dyn", "Q_dyn",
// "Qdp_dyn"): by the time ModelInitHomme::run() executes (well after
// HommeDynamics::create_requests()/set_computed_group_impl()), these Fields
// are already aliased to Homme::Context's own ElementsState/Tracers views, so
// writing into them (e.g. to set the initial condition) is visible to
// Homme's dynamical core, and updating their time stamp is visible to any
// FieldManager-registered subfield of theirs (e.g. HommeDynamics's own
// RESTART-tagged internal fields), since it is literally the same Header.
//
// Since this class has program-static lifetime (like Homme::Context), while
// the Field objects it holds must not outlive Kokkos::finalize(),
// HommeDynamics::finalize_impl() clears dyn_state_fields.
class EamxxHommeContext {
public:
  static EamxxHommeContext& singleton ();

  std::map<std::string,Field> dyn_state_fields;

private:
  EamxxHommeContext () = default;
};

} // namespace scream

#endif // EAMXX_HOMME_CONTEXT_HPP
