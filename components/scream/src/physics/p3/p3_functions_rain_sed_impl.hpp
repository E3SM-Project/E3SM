#ifndef P3_FUNCTIONS_ICE_SED_IMPL_HPP
#define P3_FUNCTIONS_ICE_SED_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 rain sedimentation functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_rain_fall_velocity(
  const Smask& qr_gt_small, const view_2d_table& vn_table, const view_2d_table& vm_table,
  const Spack& qr_incld, const Spack& rcldm, const Spack& rhofacr, Spack& nr,
  Spack& nr_incld, Spack& mu_r, Spack& lamr, Spack& V_qr, Spack& V_nr)
{
  constexpr Scalar nsmall = C::NSMALL;

  nr.set(qr_gt_small, pack::max(nr, nsmall));
  nr.set(qr_gt_small, nr_incld*rcldm);

  Table3 table;
  Spack tmp1, tmp2; //ignore
  get_rain_dsd2(qr_gt_small, qr_incld, nr, mu_r, lamr, tmp1, tmp2, rcldm);

  lookup(qr_gt_small, mu_r, lamr, table);
  // mass-weighted fall speed:
  V_qr.set(qr_gt_small,
           apply_table(qr_gt_small, vm_table, table) * rhofacr);
  // number-weighted fall speed:
  V_nr.set(qr_gt_small,
           apply_table(qr_gt_small, vn_table, table) * rhofacr);
}

} // namespace p3
} // namespace scream

#endif
