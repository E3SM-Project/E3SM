#ifndef EAMXX_WATER_TRACERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACERS_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::abs for ekat::Pack

namespace scream {
namespace water_tracers {

/*
 * Device-callable functions for water tracer ratios and ocean evaporation.
 *
 * ScalarT must be an ekat::Pack<Real,N> for some N; N=1 stands in for the
 * scalar case (the same convention as eamxx_water_isotopes_fractionation.hpp).
 * Lanes are selected with ekat::Mask, and every division or comparison that
 * could see a zero, NaN or fill value is masked first, so these functions are
 * safe in FPE-trapping builds.
 *
 * Fluxes are upward positive [kg/m2/s]: the opposite sign of the coupler's evap.
 */

// Minimum water amount for a stable ratio computation.
// From water_tracers.F90 (David Noone): smaller makes the scheme more accurate.
static constexpr Real wtrc_qmin = 1.e-22;

// The coupler fills out-of-domain ocean states with spval = 1e30
// (share/util/shr_flux_mod.F90); anything this large is not a physical value.
static constexpr Real ocean_state_fill_threshold = 1.e20;

// How a generic tracer's evaporation flux is formed (see generic_ocean_flux).
enum class OceanFlux { BulkProportional, BulkFormula };

// -----------------------------------------------------------------------
// Tracer ratio R = qtrc/qtot, or r_std where |qtot| < wtrc_qmin.
//
// Port of water_tracers.F90::wtrc_ratio (David Noone). The denominator is
// replaced by 1 on the masked lanes BEFORE dividing, so a dry cell
// (qtot = qtrc = 0) never forms 0/0.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_ratio (const ScalarT& qtrc, const ScalarT& qtot,
                    const typename ekat::ScalarTraits<ScalarT>::scalar_type r_std = 1)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto tiny = ekat::abs(qtot) < RealT(wtrc_qmin);

  ScalarT denom = qtot;
  denom.set(tiny, RealT(1));

  ScalarT ratio = qtrc / denom;
  ratio.set(tiny, r_std);
  return ratio;
}

// -----------------------------------------------------------------------
// Bulk ocean flux of a water species, upward positive:
//
//   F = f * rho * u * r * (s_eff - q_i)
//
// f = ocean fraction, rho = bottom-level air density, u = ocean friction
// velocity (So_ustar), r = moisture transfer coefficient (So_re, the square
// root of the exchange coefficient in the Large scheme), s_eff = the species'
// saturation specific humidity at the sea surface, q_i = the species at the
// bottom level. With s_eff = So_ssq and q_i = qv this is the coupler's
// evaporation (shr_flux_mod.F90: evap = rbot*ustar*re*(qbot - ssq)) with the
// opposite sign, times f because the imported states are not fraction
// weighted.
//
// Callers rely on this exact operation order: E_ocn is computed with it, so
// e.g. bulk_ocean_flux(f,rho,u,r,s*1,q) == E_ocn bit-for-bit.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT bulk_ocean_flux (const ScalarT& f, const ScalarT& rho, const ScalarT& u,
                         const ScalarT& r, const ScalarT& s_eff, const ScalarT& q_i)
{
  return f*rho*u*r*(s_eff - q_i);
}

// -----------------------------------------------------------------------
// Valid-ocean mask: f > 0, u > 0, and finite r, s with |r|, |s| < 1e20.
//
// NaN lanes are found with isnan and zeroed before any ordered comparison,
// because an ordered comparison with a NaN operand may raise FE_INVALID.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ekat::Mask<ScalarT::n> ocean_valid (const ScalarT& f, const ScalarT& u,
                                    const ScalarT& r, const ScalarT& s)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto nan = ekat::isnan(f) || ekat::isnan(u) || ekat::isnan(r) || ekat::isnan(s);

  ScalarT fz = f, uz = u, rz = r, sz = s;
  fz.set(nan, RealT(0));
  uz.set(nan, RealT(0));
  rz.set(nan, RealT(0));
  sz.set(nan, RealT(0));

  const RealT fill = ocean_state_fill_threshold;
  return !nan && fz > RealT(0) && uz > RealT(0) &&
         ekat::abs(rz) < fill && ekat::abs(sz) < fill;
}

// -----------------------------------------------------------------------
// Bulk ocean evaporation E_ocn = bulk_ocean_flux(f, rho, u, r, s, q).
//
// On return f, u, r and s are zeroed wherever the ocean state is not valid
// (see ocean_valid), so no later arithmetic touches coupler fill values, and
// E_ocn is exactly 0 there. Downstream flux functions use f > 0 as the
// valid-ocean test.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ocean_evaporation (ScalarT& f, ScalarT& u, ScalarT& r, ScalarT& s,
                           const ScalarT& rho, const ScalarT& q)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto invalid = !ocean_valid(f, u, r, s);
  f.set(invalid, RealT(0));
  u.set(invalid, RealT(0));
  r.set(invalid, RealT(0));
  s.set(invalid, RealT(0));

  ScalarT e_ocn(0);
  e_ocn.set(!invalid, bulk_ocean_flux(f, rho, u, r, s, q));
  return e_ocn;
}

// -----------------------------------------------------------------------
// Land and sea-ice tracer evaporation is out of scope for now. This hook is
// added to every mass-tracer flux so a later change only needs to fill it in.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT land_ice_tracer_flux ()
{
  return ScalarT(0);
}

// -----------------------------------------------------------------------
// Ocean flux of a generic tracer with ocean ratio r_ocn:
//
//   not valid ocean (f <= 0):  0
//   E_ocn > 0, BulkProportional: r_ocn * E_ocn
//   E_ocn > 0, BulkFormula:      bulk_ocean_flux(f, rho, u, r, s*r_ocn, q_i)
//   E_ocn <= 0 (dew):            wtrc_ratio(q_i, q) * E_ocn
//
// BulkFormula is the isotope formula with alpha = 1. Dew removes the tracer
// at its bottom-level ratio to qv.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT generic_ocean_flux (const OceanFlux kind,
                            const typename ekat::ScalarTraits<ScalarT>::scalar_type r_ocn,
                            const ScalarT& f, const ScalarT& rho, const ScalarT& u,
                            const ScalarT& r, const ScalarT& s, const ScalarT& q,
                            const ScalarT& q_i, const ScalarT& e_ocn)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto ocean = f > RealT(0);
  const auto evap  = ocean && e_ocn > RealT(0);
  const auto dew   = ocean && !(e_ocn > RealT(0));

  ScalarT flux(0);
  if (evap.any()) {
    if (kind == OceanFlux::BulkProportional) {
      flux.set(evap, r_ocn*e_ocn);
    } else {
      flux.set(evap, bulk_ocean_flux(f, rho, u, r, s*r_ocn, q_i));
    }
  }
  if (dew.any()) {
    flux.set(dew, wtrc_ratio(q_i, q)*e_ocn);
  }
  return flux;
}

// -----------------------------------------------------------------------
// Ocean flux of a property tracer p carried by the mass tracer c, where
// P = scale*source + offset and f_c is the carrier's ocean flux:
//
//   not valid ocean (f <= 0):  0
//   E_ocn > 0:                 P * f_c
//   E_ocn <= 0 (dew):          wtrc_ratio(q_p, q_c) * f_c
//
// The branch keys on the bulk water flux E_ocn, as for every other tracer
// kind, not on the sign of f_c: evaporating water carries the property of
// the source region, and dew removes the tracer at its current ratio.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT property_ocean_flux (const ScalarT& f, const ScalarT& e_ocn, const ScalarT& p,
                             const ScalarT& f_c, const ScalarT& q_p, const ScalarT& q_c)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto ocean = f > RealT(0);
  const auto evap  = ocean && e_ocn > RealT(0);
  const auto dew   = ocean && !(e_ocn > RealT(0));

  ScalarT flux(0);
  if (evap.any()) {
    flux.set(evap, p*f_c);
  }
  if (dew.any()) {
    flux.set(dew, wtrc_ratio(q_p, q_c)*f_c);
  }
  return flux;
}

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACERS_FUNCTIONS_HPP
