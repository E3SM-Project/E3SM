#ifndef ZM_ZM_CALC_FRACTIONAL_ENTRAINMENT_IMPL_HPP
#define ZM_ZM_CALC_FRACTIONAL_ENTRAINMENT_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_calc_fractional_entrainment. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_calc_fractional_entrainment(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jb, // updraft base level
  const Int& jt, // updraft top level
  // Inputs/Outputs
  Int& j0, // level where updraft begins detraining
  // Inputs
  const uview_1d<const Real>& z_mid, // env altitude at mid-point
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& dz, // layer thickness
  const uview_1d<const Real>& h_env, // env moist stat energy
  const uview_1d<const Real>& h_env_sat, // env saturated moist stat energy
  // Inputs/Outputs
  Real& h_env_min, // mid-tropospheric MSE minimum
  // Outputs
  const uview_1d<Real>& lambda, // fractional entrainment
  Real& lambda_max) // fractional entrainment maximum
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
