#ifndef SCREAM_COSP_OPTICS_IMPL_HPP
#define SCREAM_COSP_OPTICS_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

// cosp_simulator_optics2D/3D (cosp_optics.F90). The Fortran stores the subcolumn cloud
// type as a real, and compares it against the integer flags.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::simulator_optics(const S flag, const S var_cv, const S var_ls)
{
  if (flag == C::I_LSC) return var_ls;
  if (flag == C::I_CVC) return var_cv;
  return 0;
}

// MODIS_OPTICS_PARTITION (cosp_optics.F90). Pointwise; the Fortran declares the arrays as
// (npoints,nlev,ncolumns) but is called with (npoints,ncolumns,nlev) arrays, which is
// harmless since nothing couples neighbouring elements.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::modis_optics_partition(const S cloud_water, const S cloud_ice,
                                            const S water_size, const S ice_size,
                                            const S tau, S& tau_liq, S& tau_ice)
{
  S frac_liq;
  if (cloud_ice <= 0) {
    frac_liq = 1;
  } else if (cloud_water <= 0) {
    frac_liq = 0;
  } else {
    // Geometric optics limit, tau as LWP/re (proportional to LWC/re)
    frac_liq = (cloud_water/water_size) /
               (cloud_water/water_size + cloud_ice/(C::ice_density*ice_size));
  }
  tau_liq = frac_liq*tau;
  tau_ice = tau - tau_liq;
}

// modis_optics (cosp_optics.F90), for one (point, subcolumn, level). Sizes are in microns.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::modis_optics(const S tau_liq, const S size_liq,
                                  const S tau_ice, const S size_ice,
                                  S& frac_liq, S& g, S& w0)
{
  const S water_g  = get_g_nir  (C::phaseIsLiquid, size_liq);
  const S water_w0 = get_ssa_nir(C::phaseIsLiquid, size_liq);
  const S ice_g    = get_g_nir  (C::phaseIsIce,    size_ice);
  const S ice_w0   = get_ssa_nir(C::phaseIsIce,    size_ice);

  // Combine ice and water optical properties
  const S tau = tau_ice + tau_liq;
  if (tau > 0) {
    w0 = (tau_liq*water_w0 + tau_ice*ice_w0) / tau;
    g  = (tau_liq*water_g*water_w0 + tau_ice*ice_g*ice_w0) / (w0*tau);
  } else {
    w0 = 0;
    g  = 0;
  }

  // Proportion of the optical thickness due to liquid
  if (tau_liq + tau_ice > 0) {
    frac_liq = tau_liq/(tau_liq + tau_ice);
  } else {
    frac_liq = 0;
  }
}

namespace impl {

// fit_to_quadratic and fit_to_cubic (modis_simulator.F90)
template <typename S>
KOKKOS_INLINE_FUNCTION
S fit_to_quadratic(const S x, const S c1, const S c2, const S c3)
{
  return c1 + x*(c2 + x*(c3));
}

template <typename S>
KOKKOS_INLINE_FUNCTION
S fit_to_cubic(const S x, const S c1, const S c2, const S c3, const S c4)
{
  return c1 + x*(c2 + x*(c3 + x*c4));
}

} // namespace impl

// get_g_nir (modis_simulator.F90). The fit coefficients are real(wp) parameters
// initialized from single precision literals, hence the float casts.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::get_g_nir(const int phase, const S re)
{
  using impl::fit_to_quadratic;
  using impl::fit_to_cubic;

  const S ice1 = S(0.7490f), ice2 = S(6.5153e-3f), ice3 = S(-5.4136e-5f);
  const S sw1  = S(1.0364f), sw2  = S(-8.8800e-2f), sw3 = S(7.0000e-3f);
  const S bw1  = S(0.6035f), bw2  = S(2.8993e-2f), bw3  = S(-1.1051e-3f), bw4 = S(1.5134e-5f);

  S g;
  if (phase == C::phaseIsLiquid) {
    if (re < 7) {
      g = fit_to_quadratic(re, sw1, sw2, sw3);
      if (re < C::re_water_min) g = fit_to_quadratic(C::re_water_min, sw1, sw2, sw3);
    } else {
      g = fit_to_cubic(re, bw1, bw2, bw3, bw4);
      if (re > C::re_water_max) g = fit_to_cubic(C::re_water_max, bw1, bw2, bw3, bw4);
    }
  } else {
    g = fit_to_quadratic(re, ice1, ice2, ice3);
    if (re < C::re_ice_min) g = fit_to_quadratic(C::re_ice_min, ice1, ice2, ice3);
    if (re > C::re_ice_max) g = fit_to_quadratic(C::re_ice_max, ice1, ice2, ice3);
  }
  return g;
}

// get_ssa_nir (modis_simulator.F90); see get_g_nir for the float casts
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::get_ssa_nir(const int phase, const S re)
{
  using impl::fit_to_quadratic;
  using impl::fit_to_cubic;

  const S ice1 = S(0.9625f), ice2 = S(-1.8069e-2f), ice3 = S(3.3281e-4f), ice4 = S(-2.2865e-6f);
  const S w1   = S(1.0044f), w2   = S(-1.1397e-2f), w3   = S(1.3300e-4f);

  S w0;
  if (phase == C::phaseIsLiquid) {
    w0 = fit_to_quadratic(re, w1, w2, w3);
    if (re < C::re_water_min) w0 = fit_to_quadratic(C::re_water_min, w1, w2, w3);
    if (re > C::re_water_max) w0 = fit_to_quadratic(C::re_water_max, w1, w2, w3);
  } else {
    w0 = fit_to_cubic(re, ice1, ice2, ice3, ice4);
    if (re < C::re_ice_min) w0 = fit_to_cubic(C::re_ice_min, ice1, ice2, ice3, ice4);
    if (re > C::re_ice_max) w0 = fit_to_cubic(C::re_ice_max, ice1, ice2, ice3, ice4);
  }
  return w0;
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_OPTICS_IMPL_HPP
