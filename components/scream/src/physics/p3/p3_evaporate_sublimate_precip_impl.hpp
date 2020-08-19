#ifndef P3_EVAPORATE_SUBLIMATE_PRECIP_IMPL_HPP
#define P3_EVAPORATE_SUBLIMATE_PRECIP_IMPL_HPP

#include "p3_functions.hpp"
#include "physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evaporate_sublimate_precip(
  const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qi_incld,
  const Spack& cld_frac_l, const Spack& cld_frac_r, const Spack& qv_sat_l, const Spack& ab, const Spack& epsr,
  const Spack& qv, Spack& qr2qv_evap_tend, Spack& nr_evap_tend,
  const Smask& context)
{
  /* It is assumed that macrophysics handles condensation/evaporation of qc and
     that there is no condensation of rain. Thus qccon, qrcon and qcevp have
     been removed from the original P3-WRF.

     Determine temporary cloud fraction, set to zero if cloud water + ice is
     very small.  This will ensure that evap/subl of precip occurs over entire */

  constexpr Scalar QSMALL   = C::QSMALL;

  const auto summ_qc_qi_incld = qc_incld + qi_incld;
  const auto set_cld_zero = summ_qc_qi_incld < sp(1.e-6) && context;

  Spack cld;
  cld.set(set_cld_zero, 0);
  cld.set(!set_cld_zero && context, cld_frac_l);

  //Only calculate if there is some rain fraction > cloud fraction
  qr2qv_evap_tend = 0;

  const auto cld_frac_r_gt_cld = cld_frac_r > cld && context;
  const auto qr_incld_ge_qsmall = qr_incld >= QSMALL && context;

  //calculate q for out-of-cloud region
  Spack qclr;
  if(cld_frac_r_gt_cld.any()){
    qclr.set(cld_frac_r_gt_cld,(qv-cld*qv_sat_l)/(1-cld));
  }

  //rain evaporation
  if(cld_frac_r_gt_cld.any() && qr_incld_ge_qsmall.any()){
    qr2qv_evap_tend.set(cld_frac_r_gt_cld && qr_incld_ge_qsmall, epsr * (qclr-qv_sat_l)/ab);
  }

  //only evap in out-of-cloud region
  if(cld_frac_r_gt_cld.any()){
    qr2qv_evap_tend.set(cld_frac_r_gt_cld,-min(qr2qv_evap_tend*(cld_frac_r-cld),0));
    qr2qv_evap_tend.set(cld_frac_r_gt_cld,qr2qv_evap_tend/cld_frac_r);
  }

  const auto qr_incld_gt_qsmall = qr_incld > QSMALL && context;
  if(qr_incld_gt_qsmall.any()){
    nr_evap_tend.set(qr_incld_gt_qsmall, qr2qv_evap_tend*(nr_incld/qr_incld));
  }
}


} // namespace p3
} // namespace scream

#endif

