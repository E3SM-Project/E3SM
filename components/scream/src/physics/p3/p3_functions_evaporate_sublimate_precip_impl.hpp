#ifndef P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP
#define P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evaporate_sublimate_precip(const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qitot_incld,
			     const Spack& lcldm, const Spack& rcldm, const Spack& qvs, const Spack& ab, const Spack& epsr,
			     const Spack& qv, Spack& qrevp, Spack& nrevp)
{
  /* It is assumed that macrophysics handles condensation/evaporation of qc and
     that there is no condensation of rain. Thus qccon, qrcon and qcevp have
     been removed from the original P3-WRF.

     Determine temporary cloud fraction, set to zero if cloud water + ice is
     very small.  This will ensure that evap/subl of precip occurs over entire */

  constexpr Scalar QSMALL   = C::QSMALL;

  const auto summ_qc_qitot_incld = qc_incld + qitot_incld;
  const auto set_cld_zero = summ_qc_qitot_incld < sp(1.e-6);

  Spack cld;
  cld.set(set_cld_zero,0);
  cld.set(!set_cld_zero,lcldm);

  //Only calculate if there is some rain fraction > cloud fraction
  qrevp = 0;

  const auto rcldm_gt_cld = rcldm > cld;
  const auto qr_incld_ge_qsmall = qr_incld >= QSMALL;

  //calculate q for out-of-cloud region
  Spack qclr;
  qclr.set(rcldm_gt_cld,(qv-cld*qvs)/(sp(1)-cld));

  //rain evaporation
  qrevp.set(rcldm_gt_cld && qr_incld_ge_qsmall, epsr * (qclr-qvs)/ab);

  //only evap in out-of-cloud region
  qrevp.set(rcldm_gt_cld,-min(qrevp*(rcldm-cld),0));
  qrevp.set(rcldm_gt_cld,qrevp/rcldm);

  const auto qr_incld_gt_qsmall = qr_incld > QSMALL;

  nrevp.set(qr_incld_gt_qsmall, qrevp*(nr_incld/qr_incld));
}


} // namespace p3
} // namespace scream

#endif

