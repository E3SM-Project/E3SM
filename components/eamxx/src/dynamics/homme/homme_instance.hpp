#ifndef EAMXX_HOMME_INSTANCE_HPP
#define EAMXX_HOMME_INSTANCE_HPP

#include "share/field/field.hpp"
#include "share/remap/abstract_remapper.hpp"

#include <map>
#include <memory>
#include <string>

namespace scream
{

// The EAMxx-side handle on the (one) live Homme instance in this process.
//
// It owns (shares) the EAMxx Fields that are aliased by Homme's own C++ data structures (the
// dynamics-grid state, tracers, forcing, ...), and implements the steps of the init sequence
// that only involve those, with no dependency on the AtmosphereProcess interface. This way,
// different clients (e.g., the HommeDynamics atm process, and an initializer of the atm
// input fields) can operate on the same Homme data.
//
// Clients get the instance through get(). The instance lives as long as at least one client
// holds a pointer to it, and it is re-created by the next call to get() if it was destroyed.
// Since the Fields in the instance must not outlive Kokkos, clients should release their
// pointer before Kokkos is finalized.
class HommeInstance {
public:
  static std::shared_ptr<HommeInstance> get ();

  ~HommeInstance () = default;

  // Share a Field living on the dynamics grid with the instance.
  // The field must be one of those that Homme's data structures alias (see init_views).
  void add_dyn_field (const std::string& name, const Field& f);

  bool has_dyn_field (const std::string& name) const { return m_dyn_fields.count(name)>0; }

  // Fields are cheap handles to the same data, so we can return them by value
  Field get_dyn_field (const std::string& name) const;

  // Set the views of the dyn fields into the Homme c++ structures.
  // The following fields must have been added: v_dyn, vtheta_dp_dyn, phi_int_dyn, w_int_dyn,
  // dp3d_dyn, ps_dyn, omega_dyn, Q_dyn, Qdp_dyn, FQ_dyn, FT_dyn, FM_dyn. If Homme is set up
  // for 3d turbulence, also Km_dyn and Kh_dyn.
  void init_views ();

  // Copy states from the n0 timelevel to all the other timelevels
  void copy_dyn_states_to_all_timelevels ();

  // Init Homme's dynamics-grid state from the initial conditions on the GLL grid
  // (horiz_winds, T_mid, ps, phis, tracers). The pseudo density is first computed,
  // hydrostatically, on the GLL grid (so pseudo_density_gll is an output of this
  // function), and then everything is remapped to the dyn grid with the input remapper
  // (which must be a GLL->dyn remapper, with no field registered). Then, the other state
  // variables (vtheta_dp, phi_int, Qdp) are computed, and the n0 timelevel is copied to all
  // the others.
  // NOTE: this method does NOT update any time stamp.
  void init_dyn_state_from_gll_ic (AbstractRemapper& gll2dyn,
                                   const Field& horiz_winds_gll,
                                   const Field& T_mid_gll,
                                   const Field& ps_gll,
                                   const Field& phis_gll,
                                   const Field& tracers_gll,
                                   Field& pseudo_density_gll);

private:
  HommeInstance () = default;

  // The fields on the dynamics grid
  std::map<std::string,Field> m_dyn_fields;
};

} // namespace scream

#endif // EAMXX_HOMME_INSTANCE_HPP
