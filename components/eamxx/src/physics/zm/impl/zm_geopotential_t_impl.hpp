#ifndef ZM_ZM_GEOPOTENTIAL_T_IMPL_HPP
#define ZM_ZM_GEOPOTENTIAL_T_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_geopotential_t. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_geopotential_t(
  // Inputs
  const MemberType& team,
  const Int& pver,                    // number of mid-point levels
  const Int& pverp,                   // number of interface levels
  const uview_1d<const Real>& pint,   // interface pressures               [Pa]  [pverp]
  const uview_1d<const Real>& pmid,   // midpoint pressures                [Pa]  [pver]
  const uview_1d<const Real>& pdel,   // layer thickness                   [Pa]  [pver]
  const uview_1d<const Real>& t,      // temperature                       [K]   [pver]
  const uview_1d<const Real>& q,      // specific humidity                 [kg/kg] [pver]
  // Inputs/Outputs
  const uview_1d<Real>& zi,           // height above surface at interfaces [m]  [pverp]
  const uview_1d<Real>& zm)           // geopotential height at mid level   [m]  [pver]
{
  //----------------------------------------------------------------------------
  // Purpose: Compute the geopotential height (above the surface) at the
  // midpoints and interfaces using the input temperatures and pressures.
  // Copied and modified from geopotential.F90 (zm_geopotential_t).
  //----------------------------------------------------------------------------
  // Physical constants reproduced from zm_eamxx_bridge_physconst (shr_const_mod.F90).
  // These are computed here rather than taken from PC:: so that results are
  // bit-for-bit identical to the Fortran zm_geopotential_t reference, which uses
  // shr_const-derived values for rair and zvir (they differ from PC::Rair/PC::ZVIR).
  constexpr Real shr_boltz  = 1.38065e-23;          // Boltzmann's constant   [J/K/molecule]
  constexpr Real shr_avogad = 6.02214e26;           // Avogadro's number      [molecules/kmole]
  constexpr Real shr_mwdair = 28.966;               // molecular weight dry air [kg/kmole]
  constexpr Real shr_mwwv   = 18.016;               // molecular weight water vapor [kg/kmole]
  constexpr Real rgas       = shr_avogad*shr_boltz; // universal gas constant [J/K/kmole]
  constexpr Real rair       = rgas/shr_mwdair;      // dry air gas constant   [J/K/kg]
  constexpr Real rwv        = rgas/shr_mwwv;        // water vapor gas const  [J/K/kg]
  constexpr Real zvir       = (rwv/rair) - Real(1); // virtual temp. factor (rh2o/rair)-1
  constexpr Real gravit     = 9.80616;              // gravitational accel.   [m/s2]

  //----------------------------------------------------------------------------
  // The geopotential is built from the bottom up via a recurrence in k
  // (zi(k) depends on zi(k+1)), so it must run on a single thread per team.
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    // The lowest height interface is at the surface, and thus zero by definition
    zi(pverp-1) = 0;
    // Compute zi, zm from bottom up
    for (Int k = pver-1; k >= 0; --k) {
      // First set hydrostatic elements consistent with dynamics
      const Real hkl = pdel(k) / pmid(k);
      const Real hkk = 0.5 * hkl;
      // Now compute tv, zm, zi
      const Real tvfac = 1 + zvir * q(k);
      const Real tv    = t(k) * tvfac;
      zm(k) = zi(k+1) + (rair/gravit) * tv * hkk;
      zi(k) = zi(k+1) + (rair/gravit) * tv * hkl;
    }
  });
}

} // namespace zm
} // namespace scream

#endif
