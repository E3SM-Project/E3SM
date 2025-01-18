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
  const P3Runtime& runtime_options,
  const Smask& context)
{
  // ------------------------------------------------------
  // self-collection and breakup of rain
  // (breakup following modified Verlinde and Cotton scheme)

  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar rho_h2o  = C::RHO_H2O;
  constexpr Scalar pi       = C::Pi;

  const Scalar rain_selfcollection_breakup_diameter =
      runtime_options.rain_selfcollection_breakup_diameter;
  const Scalar rain_selfcollection_prefactor =
      runtime_options.rain_selfcollection_prefactor;

  const auto qr_incld_not_small = qr_incld >= qsmall && context;

  if(qr_incld_not_small.any()) {
    // use mass-mean diameter (do this by using
    // the old version of lambda w/o mu dependence)
    // note there should be a factor of 6^(1/3), but we
    // want to keep breakup threshold consistent so 'dum'
    // is expressed in terms of lambda rather than mass-mean D
    const auto dum2 = cbrt((qr_incld) / (pi * rho_h2o * nr_incld));

    Spack dum;
    const auto dum2_lt_dum1 =
        dum2 < rain_selfcollection_breakup_diameter && qr_incld_not_small;
    const auto dum2_gt_dum1 =
        dum2 >= rain_selfcollection_breakup_diameter && qr_incld_not_small;
    dum.set(dum2_lt_dum1, 1);
    if(dum2_gt_dum1.any()) {
      dum.set(dum2_gt_dum1,
              2 - exp(2300 * (dum2 - rain_selfcollection_breakup_diameter)));
    }

    nr_selfcollect_tend.set(
        qr_incld_not_small,
        dum * rain_selfcollection_prefactor * nr_incld * qr_incld * rho);
  }
}

} // namespace p3
} // namespace scream

#endif // P3_RAIN_SELF_COLLECTION_IMPL_HPP
