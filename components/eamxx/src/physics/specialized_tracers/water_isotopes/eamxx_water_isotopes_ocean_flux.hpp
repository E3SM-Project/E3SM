#ifndef EAMXX_WATER_ISOTOPES_OCEAN_FLUX_HPP
#define EAMXX_WATER_ISOTOPES_OCEAN_FLUX_HPP

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"
#include "eamxx_water_isotopes_fractionation.hpp"
#include "eamxx_water_isotopes_parameters.hpp"

namespace scream
{
namespace water_isotopes
{

/*
 * Upward (positive) ocean evaporation flux of one water isotopologue at the
 * bottom of a column [kg/m2/s], with equilibrium fractionation only:
 *
 *   f <= 0 (not a valid ocean point): F_i = 0
 *   E_ocn > 0 (evaporation):          F_i = bulk_ocean_flux(f, rho, u, r, s*R_ocn/alpha, q_i)
 *   E_ocn <= 0 (dew):                 F_i = water_tracers::tracer_ratio(q_i, q) * E_ocn
 *
 * alpha = R_liquid/R_vapor at the SST, so s*R_ocn/alpha is the isotopologue's
 * saturation specific humidity over the ocean (the delq term of iCAM's
 * wiso_flxoce). For H216O alpha == 1 exactly, and the default normalized
 * ratios give R_ocn == 1, so q_i == q yields F_i == E_ocn bit-for-bit.
 *
 * Inputs, all at the bottom level: f = ocnfrac (the caller zeroes it where the
 * ocean state is invalid), rho = air density, u = ocean friction velocity,
 * r = moisture transfer coefficient (re), s = ocean saturation specific
 * humidity (ssq), q = qv, q_i = this tracer, e_ocn = bulk_ocean_flux(f, rho,
 * u, r, s, q), t_s = SST [K], r_ocn = the species' ocean ratio.
 *
 * alpha is evaluated only on evaporating lanes. Elsewhere the SST may be
 * missing or garbage, and passing it to alpha's temperature guard would abort
 * the run. ScalarT is ekat::Pack<Real,N>; WaterIsotopes calls it with N=1.
 */
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT isotope_ocean_flux (
    const ScalarT& f, const ScalarT& rho, const ScalarT& u, const ScalarT& r,
    const ScalarT& s, const ScalarT& q, const ScalarT& q_i, const ScalarT& e_ocn,
    const ScalarT& t_s,
    const typename ekat::ScalarTraits<ScalarT>::scalar_type r_ocn,
    const WaterIsotopologues species,
    const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto ocean = f > RealT(0);
  const auto evap  = ocean && (e_ocn > RealT(0));
  const auto dew   = ocean && !(e_ocn > RealT(0));

  ScalarT flux(0);
  if (evap.any()) {
    // Lanes outside evap get an in-range placeholder temperature inside
    // alpha_liquid_vapor; their result is discarded by the masked set.
    const ScalarT alpha = WaterIsotopeFractionation::alpha_liquid_vapor(
        t_s, species, CondensedOverVapor, iso_params, evap);
    flux.set(evap, water_tracers::bulk_ocean_flux(f, rho, u, r, s*r_ocn/alpha, q_i));
  }
  if (dew.any()) {
    flux.set(dew, water_tracers::tracer_ratio(q_i, q)*e_ocn);
  }
  return flux;
}

} // namespace water_isotopes
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_OCEAN_FLUX_HPP
