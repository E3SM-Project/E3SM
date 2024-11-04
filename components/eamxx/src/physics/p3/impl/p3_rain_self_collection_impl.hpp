#ifndef P3_RAIN_SELF_COLLECTION_IMPL_HPP
#define P3_RAIN_SELF_COLLECTION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_self_collection(
  const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nr_selfcollect_tend,
  const physics::P3_Constants<S> & p3constants,
  const Smask& context)
{
  // ------------------------------------------------------
  // self-collection and breakup of rain
  // (breakup following modified Verlinde and Cotton scheme)

  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar rho_h2o  = C::RHO_H2O;
  constexpr Scalar pi       = C::Pi;

  const Scalar p3_d_breakup_cutoff = p3constants.p3_d_breakup_cutoff;

  const auto qr_incld_not_small = qr_incld >= qsmall && context;

  if (qr_incld_not_small.any()) {
    
    const auto dum2 = cbrt((qr_incld)/(pi*rho_h2o*nr_incld));

    Spack dum;
    const auto dum2_lt_dum1 = dum2 < p3_d_breakup_cutoff && qr_incld_not_small;
    const auto dum2_gt_dum1 = dum2 >= p3_d_breakup_cutoff && qr_incld_not_small;
    dum.set(dum2_lt_dum1, 1);
    if (dum2_gt_dum1.any()) {
      dum.set(dum2_gt_dum1, 2 - exp(2300 * (dum2-p3_d_breakup_cutoff)));
    }

    nr_selfcollect_tend.set(qr_incld_not_small, dum*sp(5.78)*nr_incld*qr_incld*rho);
  }
}

} // namespace p3
} // namespace scream

#endif // P3_RAIN_SELF_COLLECTION_IMPL_HPP
