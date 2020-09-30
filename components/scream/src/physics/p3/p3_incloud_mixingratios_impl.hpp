#ifndef P3_INCLOUD_MIXINGRATIOS_IMPL_HPP
#define P3_INCLOUD_MIXINGRATIOS_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calculate_incloud_mixingratios(
  const Spack& qc, const Spack& qr, const Spack& qi, const Spack& qm, const Spack& nc,
  const Spack& nr, const Spack& ni, const Spack& bm, const Spack& inv_cld_frac_l,
  const Spack& inv_cld_frac_i, const Spack& inv_cld_frac_r,
  Spack& qc_incld, Spack& qr_incld, Spack& qi_incld, Spack& qm_incld,
  Spack& nc_incld, Spack& nr_incld, Spack& ni_incld, Spack& bm_incld,
  const Smask& context)
{
   constexpr Scalar qsmall = C::QSMALL;
   constexpr Scalar incloud_limit = C::incloud_limit;
   constexpr Scalar precip_limit  = C::precip_limit;

   const auto qc_ge_qsmall = qc >= qsmall && context;
   const auto not_qc_ge_qsmall = !qc_ge_qsmall && context;

   const auto qi_ge_qsmall = qi >= qsmall && context;
   const auto not_qi_ge_qsmall = !qi_ge_qsmall && context;

   const auto qi_and_qm_ge_qsmall = (qi >= qsmall) && (qm >= qsmall) && context;
   const auto not_qi_and_qm_ge_qsmall = !qi_and_qm_ge_qsmall  && context;

   const auto qr_ge_qsmall = qr >= qsmall && context;
   const auto not_qr_ge_qsmall = !qr_ge_qsmall && context;

   qc_incld.set(qc_ge_qsmall, qc*inv_cld_frac_l);
   nc_incld.set(qc_ge_qsmall, max(nc*inv_cld_frac_l, 0));

   qc_incld.set(not_qc_ge_qsmall, 0);
   nc_incld.set(not_qc_ge_qsmall, 0);

   qi_incld.set(qi_ge_qsmall, qi*inv_cld_frac_i);
   ni_incld.set(qi_ge_qsmall, max(ni*inv_cld_frac_i, 0));

   qi_incld.set(not_qi_ge_qsmall, 0);
   ni_incld.set(not_qi_ge_qsmall, 0);

   qm_incld.set(qi_and_qm_ge_qsmall, qm*inv_cld_frac_i);
   bm_incld.set(qi_and_qm_ge_qsmall, max(bm*inv_cld_frac_l, 0));

   qm_incld.set(not_qi_and_qm_ge_qsmall, 0);
   bm_incld.set(not_qi_and_qm_ge_qsmall, 0);

   qr_incld.set(qr_ge_qsmall, qr*inv_cld_frac_r);
   nr_incld.set(qr_ge_qsmall, max(nr*inv_cld_frac_r, 0));

   qr_incld.set(not_qr_ge_qsmall, 0);
   nr_incld.set(not_qr_ge_qsmall, 0);

   const auto any_gt_limit =
     (qc_incld > incloud_limit || qi_incld > incloud_limit || qr_incld > precip_limit || bm_incld > incloud_limit) && context;

   qc_incld.set(any_gt_limit, min(qc_incld, incloud_limit));
   qi_incld.set(any_gt_limit, min(qi_incld, incloud_limit));
   bm_incld.set(any_gt_limit, min(bm_incld, incloud_limit));
   qr_incld.set(any_gt_limit, min(qr_incld, precip_limit));
}

} // namespace p3
} // namespace scream

#endif
