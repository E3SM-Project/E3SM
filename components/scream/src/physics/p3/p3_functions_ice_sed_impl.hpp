#ifndef P3_FUNCTIONS_ICE_SED_IMPL_HPP
#define P3_FUNCTIONS_ICE_SED_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>
::calc_bulk_rho_rime(
  const Smask& qi_gt_small, const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim)
{
  constexpr Scalar bsmall       = C::BSMALL;
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar rho_rime_min = C::rho_rimeMin;
  constexpr Scalar rho_rime_max = C::rho_rimeMax;

  Spack rho_rime(0);

  Smask bi_rim_gt_small = bi_rim >= bsmall;
  rho_rime.set(qi_gt_small && bi_rim_gt_small, qi_rim / bi_rim);

  Smask rho_rime_lt_min = rho_rime < rho_rime_min;
  Smask rho_rime_gt_max = rho_rime > rho_rime_max;

  // impose limits on rho_rime;  adjust bi_rim if needed
  rho_rime.set(qi_gt_small && bi_rim_gt_small && rho_rime_lt_min, rho_rime_min);
  rho_rime.set(qi_gt_small && bi_rim_gt_small && rho_rime_gt_max, rho_rime_max);
  bi_rim.set(qi_gt_small && bi_rim_gt_small && (rho_rime_gt_max || rho_rime_lt_min), qi_rim / rho_rime);

  qi_rim.set(qi_gt_small && !bi_rim_gt_small, 0);
  bi_rim.set(qi_gt_small && !bi_rim_gt_small, 0);
  rho_rime.set(qi_gt_small && !bi_rim_gt_small, 0);

  // set upper constraint qi_rim <= qi_tot
  Smask qi_rim_gt_qi = (qi_rim > qi_tot) && (rho_rime > 0) && qi_gt_small;
  qi_rim.set(qi_rim_gt_qi, qi_tot);
  bi_rim.set(qi_rim_gt_qi, qi_rim / rho_rime);

  // impose consistency
  Smask qi_rim_lt_small = (qi_rim < qsmall) && qi_gt_small;
  qi_rim.set(qi_rim_lt_small, 0);
  bi_rim.set(qi_rim_lt_small, 0);

  return rho_rime;
}

} // namespace p3
} // namespace scream

#endif
