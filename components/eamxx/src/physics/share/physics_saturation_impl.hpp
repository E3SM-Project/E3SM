#ifndef PHYSICS_SATURATION_IMPL_HPP
#define PHYSICS_SATURATION_IMPL_HPP

#include "physics_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace physics {

/*
 * Implementation of saturation functions. Clients should NOT #include
 * this file, #include physics_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void  Functions<S,D>
::check_temperature(const Spack& t_atm, const char* caller, const Smask& range_mask)
{

  //"range_mask" masks out packs which are padded with uninitialized values

  /*Future work:
  Make error message dynamic so that it tells the user the function name (func_name) from where
  this error is coming from. Currently func_name is not used in this function */

  //NOTE: EKAT_KERNEL_REQUIRE_MSG requires first argument to be "False" to exit with an error message

  //find out if there are any negative temperatures in the pack
  const auto is_neg_t_atm = (t_atm <= 0) && range_mask;
  EKAT_KERNEL_REQUIRE_MSG(!(is_neg_t_atm.any()), caller); //exit with an error message

  //find out if there are any NaN temperatures in the pack
  const auto is_nan_t_atm = isnan(t_atm) && range_mask;
  EKAT_KERNEL_REQUIRE_MSG(!(is_nan_t_atm.any()), caller); //exit with an error message

}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::MurphyKoop_svp(const Spack& t_atm, const bool ice, const Smask& range_mask, const char* caller)
{

  //First check if the temperature is legitimate or not
  check_temperature(t_atm, caller ? caller : "MurphyKoop_svp", range_mask);

  //Formulas used below are from the following paper:
  //Murphy, D. M., and T. Koop (2005), Review of the vapour pressures of ice
  //and supercooled water for atmospheric applications, Quart J. Roy. Meteor.
  //Soc., 131(608), 1539–1565,

  Spack result;
  static constexpr  auto tmelt = C::Tmelt;
  const Smask ice_mask = (t_atm < tmelt) && ice;
  const Smask liq_mask = !ice_mask;

  if (ice_mask.any()) {
    //Equation (7) of the paper
    // (good down to 110 K)
    //creating array for storing coefficients of ice sat equation
    static constexpr Scalar ic[]= {9.550426, 5723.265, 3.53068, 0.00728332};
    const Spack ice_result = exp(ic[0] - (ic[1] / t_atm) + (ic[2] * log(t_atm)) - (ic[3] * t_atm));

    result.set(ice_mask, ice_result);
  }

  if (liq_mask.any()) {
    //Equation (10) of the paper
    // (good for 123 < T < 332 K)
    //creating array for storing coefficients of liq sat equation
    static constexpr Scalar lq[] = {54.842763, 6763.22, 4.210, 0.000367, 0.0415, 218.8, 53.878,
			 1331.22, 9.44523, 0.014025 };
    const auto logt = log(t_atm);
    const Spack liq_result = exp(lq[0] - (lq[1] / t_atm) - (lq[2] * logt) + (lq[3] * t_atm) +
			   (tanh(lq[4] * (t_atm - lq[5])) * (lq[6] - (lq[7] / t_atm) -
							 (lq[8] * logt) + lq[9] * t_atm)));

    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::polysvp1(const Spack& t, const bool ice, const Smask& range_mask, const char* caller)
{
  // REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

  //First check if the temperature is legitimate or not
  check_temperature(t, caller ? caller : "polysvp1", range_mask);

  const Spack dt = max(t - sp(273.15), sp(-80));
  Spack result;
  static constexpr  auto tmelt = C::Tmelt;
  const Smask ice_mask = (t < tmelt) && ice;
  const Smask liq_mask = !ice_mask;

  // -------------------------------------------

  // Flatau formulation:
  if (ice_mask.any()) {
    // ice
    static constexpr Scalar ai[] = {
      6.11147274,     0.503160820,     0.188439774e-1,
      0.420895665e-3, 0.615021634e-5,  0.602588177e-7,
      0.385852041e-9, 0.146898966e-11, 0.252751365e-14};
    const Spack ice_result = (ai[0] + dt*(ai[1]+dt*(ai[2]+dt*(ai[3]+dt*(ai[4]+dt*(ai[5]+dt*(ai[6]+dt*(ai[7]+
                                                                                                ai[8]*dt))))))))*100;
    result.set(ice_mask, ice_result);
  }
  if (liq_mask.any()) {
    // liquid, V1.7
    static constexpr Scalar a[] = {
      6.11239921,      0.443987641,     0.142986287e-1,
      0.264847430e-3,  0.302950461e-5,  0.206739458e-7,
      0.640689451e-10,-0.952447341e-13,-0.976195544e-15};

    const Spack liq_result = (a[0] + dt*(a[1]+dt*(a[2]+dt*(a[3]+dt*(a[4]+dt*(a[5]+dt*(a[6]+dt*(a[7]+a[8]*dt))))))))*100;
    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice, const Smask& range_mask, const SaturationFcn func_idx, const char* caller)
{
  /*Arguments:
    ----------
  t_atm: temperature; p_atm: pressure; ice: logical for ice

  range_mask: is a mask which masks out padded values in the packs, which are uninitialized
  func_idx is an optional argument to decide which scheme is to be called for saturation vapor pressure
  Currently default is set to "MurphyKoop_svp"
  func_idx = Polysvp1 (=0) --> polysvp1 (Flatau et al. 1992)
  func_idx = MurphyKoop (=1) --> MurphyKoop_svp (Murphy, D. M., and T. Koop 2005)*/

  Spack e_pres; // saturation vapor pressure [Pa]

  switch (func_idx){
    case Polysvp1:
      e_pres = polysvp1(t_atm, ice, range_mask, caller);
      break;
    case MurphyKoop:
      e_pres = MurphyKoop_svp(t_atm, ice, range_mask, caller);
      break;
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Invalid func_idx supplied to qv_sat.");
    }

  static constexpr  auto ep_2 = C::ep_2;
  return ep_2 * e_pres / max(p_atm-e_pres, sp(1.e-3));
}

} // namespace physics
} // namespace scream

#endif // PHYSICS_SATURATION_IMPL_HPP
