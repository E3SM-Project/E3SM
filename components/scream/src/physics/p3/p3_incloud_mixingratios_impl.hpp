#ifndef P3_INCLOUD_MIXINGRATIOS_IMPL_HPP
#define P3_INCLOUD_MIXINGRATIOS_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calculate_incloud_mixingratios(
  const Spack& qc, const Spack& qr, const Spack& qitot, const Spack& qirim, const Spack& nc,
  const Spack& nr, const Spack& nitot, const Spack& birim, const Spack& inv_lcldm,
  const Spack& inv_icldm, const Spack& inv_rcldm,
  Spack& qc_incld, Spack& qr_incld, Spack& qitot_incld, Spack& qirim_incld,
  Spack& nc_incld, Spack& nr_incld, Spack& nitot_incld, Spack& birim_incld,
  const Smask& context)
{
   constexpr Scalar qsmall = C::QSMALL;
   constexpr Scalar incloud_limit = C::incloud_limit;
   constexpr Scalar precip_limit  = C::precip_limit;

   const auto qc_ge_qsmall = qc >= qsmall && context;
   const auto not_qc_ge_qsmall = !qc_ge_qsmall && context;

   const auto qitot_ge_qsmall = qitot >= qsmall && context;
   const auto not_qitot_ge_qsmall = !qitot_ge_qsmall && context;

   const auto qitot_and_qirim_ge_qsmall = (qitot >= qsmall) && (qirim >= qsmall) && context;
   const auto not_qitot_and_qirim_ge_qsmall = !qitot_and_qirim_ge_qsmall  && context;

   const auto qr_ge_qsmall = qr >= qsmall && context;
   const auto not_qr_ge_qsmall = !qr_ge_qsmall && context;

   qc_incld.set(qc_ge_qsmall, qc*inv_lcldm);
   nc_incld.set(qc_ge_qsmall, pack::max(nc*inv_lcldm, 0));

   qc_incld.set(not_qc_ge_qsmall, 0);
   nc_incld.set(not_qc_ge_qsmall, 0);

   qitot_incld.set(qitot_ge_qsmall, qitot*inv_icldm);
   nitot_incld.set(qitot_ge_qsmall, pack::max(nitot*inv_icldm, 0));

   qitot_incld.set(not_qitot_ge_qsmall, 0);
   nitot_incld.set(not_qitot_ge_qsmall, 0);

   qirim_incld.set(qitot_and_qirim_ge_qsmall, qirim*inv_icldm);
   birim_incld.set(qitot_and_qirim_ge_qsmall, pack::max(birim*inv_lcldm, 0));

   qirim_incld.set(not_qitot_and_qirim_ge_qsmall, 0);
   birim_incld.set(not_qitot_and_qirim_ge_qsmall, 0);

   qr_incld.set(qr_ge_qsmall, qr*inv_rcldm);
   nr_incld.set(qr_ge_qsmall, pack::max(nr*inv_rcldm, 0));

   qr_incld.set(not_qr_ge_qsmall, 0);
   nr_incld.set(not_qr_ge_qsmall, 0);

   const auto any_gt_limit =
     (qc_incld > incloud_limit || qitot_incld > incloud_limit || qr_incld > precip_limit || birim_incld > incloud_limit) && context;

   qc_incld.set(any_gt_limit,
                pack::min(qc_incld, incloud_limit));
   qitot_incld.set(any_gt_limit,
                   pack::min(qitot_incld, incloud_limit));
   birim_incld.set(any_gt_limit,
                   pack::min(birim_incld, incloud_limit));
   qr_incld.set(any_gt_limit,
                pack::min(qr_incld, precip_limit));
}

} // namespace p3
} // namespace scream

#endif
