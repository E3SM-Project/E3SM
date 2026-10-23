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

  // Finalize the restart of Homme's dynamics-grid state, assuming the n0 time level of the
  // states (and the n0_qdp one of Qdp) was read from the restart file: copy it to all the
  // other time levels, and recompute Q (the end-of-step tracers mixing ratio) as Qdp/dp
  // (only Qdp is in the restart file, to save space).
  void finish_dyn_state_restart ();

  // Setup the remapper between the dyn (GLL) grid and the FV physics grid, with pgN x pgN
  // physics cells per element. Must be called after Homme's init1 phase has completed.
  void init_fv_phys (const int pgN);

  // The N in pgN, or -1 if the FV physics grid is not in use
  int get_fv_phys_pgN () const { return m_pgN; }

  // Remap the Homme state at the current n0 time level to the FV physics grid.
  // All the input fields are outputs of this method, and they live on the FV physics grid.
  //  - ps, phis: (COL)
  //  - T, omega, dp: (COL,LEV)
  //  - uv: (COL,CMP=2,LEV)
  //  - q: (COL,CMP,LEV), with one component per tracer
  //  - strain3d_fv: (COL,CMP=6,LEV). Only used (and must be valid) if 3d turbulence is on.
  // The method also needs the field shear_strain3d_components_dyn, if 3d turbulence is on.
  void remap_dyn_to_fv_phys (const Field& ps, const Field& phis,
                             const Field& T, const Field& omega,
                             const Field& uv, const Field& q,
                             const Field& dp, const Field& strain3d_fv = Field());

  // Remap the tendency/state fields from the FV physics grid to the dyn grid, into Homme's
  // forcing data structures, and DSS the result.
  //  - FT: (COL,LEV), the temperature to use to compute T tendencies
  //  - FM: (COL,CMP=2,LEV), the horiz winds to use to compute the momentum tendencies
  //  - q: (COL,CMP,LEV), with one component per tracer
  //  - Km, Kh: (COL,LEV), the eddy diffusivities. Only used if 3d turbulence is on.
  void remap_fv_phys_to_dyn (const Field& FT, const Field& FM, const Field& q,
                             const Field& Km = Field(), const Field& Kh = Field());

private:
  HommeInstance () = default;

  // The N in pgN. Negative if the FV physics grid is not in use.
  int m_pgN = -1;

  // The fields on the dynamics grid
  std::map<std::string,Field> m_dyn_fields;
};

} // namespace scream

#endif // EAMXX_HOMME_INSTANCE_HPP
