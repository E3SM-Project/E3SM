#ifndef ZM_ZM_STATE_UPDATE_IMPL_HPP
#define ZM_ZM_STATE_UPDATE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_state_update. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_state_update(
  // Inputs
  const MemberType& team,
  const Int& pver,                     // number of mid-point levels
  const Int& pverp,                    // number of interface levels
  const Real& dt,                      // time step                        [s]
  const uview_1d<const Real>& pmid,    // mid-point pressure               [Pa]    [pver]
  const uview_1d<const Real>& pint,    // interface pressure               [Pa]    [pverp]
  const uview_1d<const Real>& pdel,    // pressure thickness               [Pa]    [pver]
  const uview_1d<const Real>& ptend_s, // tendency of dry static energy    [J/kg/s][pver]
  const uview_1d<const Real>& ptend_q, // tendency of water vapor          [kg/kg/s][pver]
  // Inputs/Outputs
  const uview_1d<Real>& zm,            // altitude at mid-levels           [m]     [pver]
  const uview_1d<Real>& zi,            // altitude at interfaces           [m]     [pverp]
  const uview_1d<Real>& t,             // temperature                      [K]     [pver]
  const uview_1d<Real>& qv)            // water vapor                      [kg/kg] [pver]
{
  //----------------------------------------------------------------------------
  // Purpose: Apply ZM convective tendencies to the state and recompute the
  // geopotential height. This combines the functionality of physics_update()
  // [physics_update_mod.F90] and physics_update_main() [physics_types.F90].
  //
  // Note: the surface geopotential (state_phis) and the dry static energy
  // update are intentionally omitted, matching the active code path used for
  // EAMxx (the DSE update is disabled there). This mirrors the Fortran
  // zm_state_update in zm_eamxx_bridge_methods.F90.
  //----------------------------------------------------------------------------
  // specific heat of dry air [J/K/kg] - matches zm_eamxx_bridge_physconst (shr_const_cpdair)
  constexpr Real cpair = 1004.64;

  //----------------------------------------------------------------------------
  // Apply tendencies to the local state (level-independent)
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&] (const Int& k) {
    // update water vapor
    qv(k) = qv(k) + ptend_q(k) * dt;
    // update temperature - assume that dS is really dEn, En=enthalpy=c_p*T,
    // then dT = dEn/c_p, so state%t += ds/c_p.
    t(k) = t(k) + ptend_s(k)/cpair * dt;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Recompute geopotential height from the updated temperature and humidity
  zm_geopotential_t(team, pver, pverp, pint, pmid, pdel, t, qv, zi, zm);
}

} // namespace zm
} // namespace scream

#endif
