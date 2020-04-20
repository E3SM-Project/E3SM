#ifndef P3_FUNCTIONS_MATH_IMPL_HPP
#define P3_FUNCTIONS_MATH_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::polysvp1(const Spack& t, const bool ice)
{
  // REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

  // ice
  static Scalar ai[] = {
    6.11147274,     0.503160820,     0.188439774e-1,
    0.420895665e-3, 0.615021634e-5,  0.602588177e-7,
    0.385852041e-9, 0.146898966e-11, 0.252751365e-14};

  // liquid, V1.7
  static Scalar a[] = {
    6.11239921,      0.443987641,     0.142986287e-1,
    0.264847430e-3,  0.302950461e-5,  0.206739458e-7,
    0.640689451e-10,-0.952447341e-13,-0.976195544e-15};

  Spack dt = pack::max(t - sp(273.15), sp(-80.0));
  Spack result;
  const auto tmelt = C::Tmelt;
  Smask ice_mask = (t < tmelt) && ice;
  Smask liq_mask = !ice_mask;

  // -------------------------------------------

  // Flatau formulation:
  if (ice_mask.any()) {
    Spack ice_result = (ai[0] + dt*(ai[1]+dt*(ai[2]+dt*(ai[3]+dt*(ai[4]+dt*(ai[5]+dt*(ai[6]+dt*(ai[7]+
                                                                                                ai[8]*dt))))))))*100;
    result.set(ice_mask, ice_result);
  }
  if (liq_mask.any()) {
    Spack liq_result = (a[0] + dt*(a[1]+dt*(a[2]+dt*(a[3]+dt*(a[4]+dt*(a[5]+dt*(a[6]+dt*(a[7]+a[8]*dt))))))))*100;
    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice)
{
  Spack e_pres; // saturation vapor pressure [Pa]

  e_pres = polysvp1(t_atm, ice);
  const auto ep_2 = C::ep_2;
  return ep_2 * e_pres / pack::max(p_atm-e_pres, sp(1.e-3));
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::qvsat_exact(const Spack& t_atm, const Spack& p_atm, const bool ice)
{
  // Compute saturation mixing ratio qs analytically. 
  // Derivation: Clausius Clapeyron says dqs/dT = L*qs/(Rv*T**2.).
  // Taking the integral between some T0 and a particular T yields:
  // qs(T) = qs(T0)*exp{-L/Rv*(1/T - 1/T0)}.
  // Note that pressure dependency comes from qs(T0), which can
  // be computed by using an experimentally-determined value for saturation 
  // vapor pressure (es_T0) at 273.15 K and converting to qs using
  // the exact formula qs(p) = es_T0/(pres - es_T0).

  const auto tmelt = C::Tmelt;
  const auto RV = C::RV;
  const auto LatVap = C::LatVap;
  const auto LatIce = C::LatIce;
  const auto ep_2 = C::ep_2;
  
  Spack result;
  Spack es_T0; //saturation vapor pressure at T=tmelt
  
  Smask ice_mask = (t_atm < tmelt) && ice;
  Smask liq_mask = !ice_mask;

  // -------------------------------------------
  // Note on implementation: would be more efficient to create Spacks for es_T0 and 
  // L which vary based on whether to calculate with respect to liquid or ice, then
  // to perform the conversion to qs_T0 and to get result 1x?

  if (ice_mask.any()) {
    Scalar es_T0=sp(611.147274); //taken from Flatau polynomial fit... could do better.
    Spack qs_T0 = ep_2 * es_T0 / pack::max(p_atm-es_T0, sp(1.e-3));
    Spack ice_result= qs_T0*pack::exp( -(LatVap+LatIce)/RV*(1/t_atm - 1/tmelt) );
    result.set(ice_mask, ice_result);
  }
  if (liq_mask.any()) {
    Scalar es_T0=sp(611.239921); //taken from Flatau polynomial fit... could do better.
    Spack qs_T0 = ep_2 * es_T0 / pack::max(p_atm-es_T0, sp(1.e-3));
    Spack liq_result= qs_T0*pack::exp( -LatVap/RV*(1/t_atm - 1/tmelt) );
    result.set(liq_mask, liq_result);
  }

  return result;
} // end qvsat_exact




} // namespace p3
} // namespace scream

#endif
