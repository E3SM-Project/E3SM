#ifndef P3_UPDATE_PROGNOSTICS_IMPL_HPP
#define P3_UPDATE_PROGNOSTICS_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_prognostic_ice(
  const Spack& qc2qi_hetero_freeze_tend, const Spack& qc2qi_collect_tend,  const Spack& qc2qr_ice_shed_tend, const Spack& nc_collect_tend,
  const Spack& nc2ni_immers_freeze_tend, const Spack& ncshdc, const Spack& qr2qi_collect_tend, const Spack& nr_collect_tend,
  const Spack& qr2qi_immers_freeze_tend, const Spack& nr2ni_immers_freeze_tend, const Spack& nr_ice_shed_tend, const Spack& qi2qr_melt_tend,
  const Spack& ni2nr_melt_tend, const Spack& qi2qv_sublim_tend, const Spack& qv2qi_vapdep_tend, const Spack& qv2qi_nucleat_tend,
  const Spack& ni_nucleat_tend, const Spack& ni_selfcollect_tend, const Spack& ni_sublim_tend, const Spack& qc2qi_berg_tend,
  const Spack& inv_exner, const bool do_predict_nc,
  const Smask& log_wetgrowth, const Scalar dt,  const Scalar& nmltratio, const Spack& rho_qm_cloud,
  Spack& ncheti_cnt, Spack& nicnt, Spack& ninuc_cnt, Spack& qcheti_cnt, Spack& qicnt, Spack& qinuc_cnt,
  Spack& th_atm, Spack& qv, Spack& qi, Spack& ni, Spack& qm, Spack& bm, Spack& qc,
  Spack& nc, Spack& qr, Spack& nr, const bool& use_hetfrz_classnuc,
  const Smask& context)
{
  constexpr Scalar QSMALL          = C::QSMALL;
  constexpr Scalar INV_RHO_RIMEMAX = C::INV_RHO_RIMEMAX;
  constexpr Scalar latvap          = C::LatVap;
  constexpr Scalar latice          = C::LatIce;

  if(use_hetfrz_classnuc){
    qc.set(context, qc + (-qcheti_cnt-qicnt-qc2qi_collect_tend-qc2qr_ice_shed_tend-qc2qi_berg_tend)*dt);
  }
  else{
    qc.set(context, qc + (-qc2qi_hetero_freeze_tend-qc2qi_collect_tend-qc2qr_ice_shed_tend-qc2qi_berg_tend)*dt);
  }


  if ( do_predict_nc ){
    if(use_hetfrz_classnuc){
      nc.set(context, nc + (-nc_collect_tend-ncheti_cnt-nicnt)*dt);
    }
    else{
      nc.set(context, nc + (-nc_collect_tend-nc2ni_immers_freeze_tend)*dt);
    }
  }

  qr.set(context, qr + (-qr2qi_collect_tend+qi2qr_melt_tend-qr2qi_immers_freeze_tend+qc2qr_ice_shed_tend)*dt);

  //apply factor to source for rain number from melting of ice, (ad-hoc
  // but accounts for rapid evaporation of small melting ice particles)
  nr.set(context, nr + (-nr_collect_tend-nr2ni_immers_freeze_tend+nmltratio*ni2nr_melt_tend+nr_ice_shed_tend+ncshdc)*dt);

  const auto qi_not_small = qi >= QSMALL && context;

  if ( qi_not_small.any() ) {
    bm.set(qi_not_small, bm - ((qi2qv_sublim_tend + qi2qr_melt_tend) / qi) * dt * bm);
    qm.set(qi_not_small, qm - ((qi2qv_sublim_tend + qi2qr_melt_tend) * qm / qi) * dt);
    qi.set(qi_not_small, qi - (qi2qv_sublim_tend + qi2qr_melt_tend) * dt);
  }

  if(use_hetfrz_classnuc){
    const auto dum = (qr2qi_collect_tend + qc2qi_collect_tend + qr2qi_immers_freeze_tend + qcheti_cnt+qicnt) * dt;
    qi.set(context, qi + (qv2qi_vapdep_tend + qv2qi_nucleat_tend + qc2qi_berg_tend+qinuc_cnt)*dt + dum);
    qm.set(context, qm + dum);
    bm.set(context, bm + (qr2qi_collect_tend * INV_RHO_RIMEMAX + qc2qi_collect_tend / rho_qm_cloud + (qr2qi_immers_freeze_tend +
                                                                              qcheti_cnt+qicnt) * INV_RHO_RIMEMAX) * dt);
    ni.set(context, ni + (ni_nucleat_tend - ni2nr_melt_tend - ni_sublim_tend - ni_selfcollect_tend + nr2ni_immers_freeze_tend +ncheti_cnt+nicnt+ninuc_cnt)*dt);
  }
  else{
    const auto dum = (qr2qi_collect_tend + qc2qi_collect_tend + qr2qi_immers_freeze_tend + qc2qi_hetero_freeze_tend) * dt;
    qi.set(context, qi + (qv2qi_vapdep_tend + qv2qi_nucleat_tend + qc2qi_berg_tend) * dt + dum);
    qm.set(context, qm + dum);
    bm.set(context, bm + (qr2qi_collect_tend * INV_RHO_RIMEMAX + qc2qi_collect_tend / rho_qm_cloud + (qr2qi_immers_freeze_tend +
                                                                              qc2qi_hetero_freeze_tend) * INV_RHO_RIMEMAX) * dt);
    ni.set(context, ni + (ni_nucleat_tend - ni2nr_melt_tend - ni_sublim_tend - ni_selfcollect_tend + nr2ni_immers_freeze_tend + nc2ni_immers_freeze_tend) * dt);
  }

  //PMC nCat deleted interactions_loop

  const auto qm_lt_thresh = qm < 0 && context;
  if (qm_lt_thresh.any()){
    qm.set(qm_lt_thresh, 0);
    bm.set(qm_lt_thresh, 0);
  }

  // densify under wet growth
  // -- to be removed post-v2.1.  Densification automatically happens
  //    during wet growth due to parameterized rime density --

  qm.set(log_wetgrowth && context, qi);
  bm.set(log_wetgrowth && context, qm * INV_RHO_RIMEMAX);

  // densify in above freezing conditions and melting
  // -- future work --
  //   Ideally, this will be treated with the predicted liquid fraction in ice.
  //   Alternatively, it can be simplified by tending qm -- qi
  //   and bm such that rho_rim (qm/bm) --> rho_liq during melting.
  // ==

  constexpr Scalar INV_CP = C::INV_CP;
  if(use_hetfrz_classnuc){
    qv.set(context, qv + (-qv2qi_vapdep_tend+qi2qv_sublim_tend-qv2qi_nucleat_tend-qinuc_cnt)*dt);
    th_atm.set(context, th_atm + inv_exner * ((qv2qi_vapdep_tend - qi2qv_sublim_tend + qv2qi_nucleat_tend+qinuc_cnt) * (latvap+latice) * INV_CP +
                                (qr2qi_collect_tend + qc2qi_collect_tend + qcheti_cnt+qicnt + qr2qi_immers_freeze_tend -
                                qi2qr_melt_tend + qc2qi_berg_tend) * latice * INV_CP) * dt);
  }
  else{
    qv.set(context, qv + (-qv2qi_vapdep_tend+qi2qv_sublim_tend-qv2qi_nucleat_tend)*dt);
    th_atm.set(context, th_atm + inv_exner * ((qv2qi_vapdep_tend - qi2qv_sublim_tend + qv2qi_nucleat_tend) * (latvap+latice) * INV_CP +
                                (qr2qi_collect_tend + qc2qi_collect_tend + qc2qi_hetero_freeze_tend + qr2qi_immers_freeze_tend -
                                qi2qr_melt_tend + qc2qi_berg_tend) * latice * INV_CP) * dt);
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_prognostic_liquid(
  const Spack& qc2qr_accret_tend, const Spack& nc_accret_tend,
  const Spack& qc2qr_autoconv_tend,const Spack& nc2nr_autoconv_tend, const Spack& ncautr,
  const Spack& nc_selfcollect_tend, const Spack& qr2qv_evap_tend, const Spack& nr_evap_tend, const Spack& nr_selfcollect_tend,
  const bool do_predict_nc, const bool do_prescribed_CCN, const Spack& inv_rho, const Spack& inv_exner,
  const Scalar dt, Spack& th_atm, Spack& qv, Spack& qc, Spack& nc, Spack& qr, Spack& nr,
  const Smask& context)
{
  constexpr Scalar NCCNST = C::NCCNST;
  constexpr int IPARAM    = C::IPARAM;
  constexpr Scalar INV_CP = C::INV_CP;
  constexpr Scalar latvap       = C::LatVap;

  qc.set(context, qc + (-qc2qr_accret_tend-qc2qr_autoconv_tend)*dt);
  qr.set(context, qr + (qc2qr_accret_tend+qc2qr_autoconv_tend-qr2qv_evap_tend)*dt);

  if (do_predict_nc || do_prescribed_CCN) {
    nc.set(context, nc + (-nc_accret_tend-nc2nr_autoconv_tend+nc_selfcollect_tend)*dt);
  }
  else {
    nc.set(context, NCCNST * inv_rho);
  }

  if (IPARAM == 1 || IPARAM == 2) {
    nr.set(context, nr + (sp(0.5) * nc2nr_autoconv_tend - nr_selfcollect_tend - nr_evap_tend) * dt);
  }
  else {
    nr.set(context, nr + (ncautr - nr_selfcollect_tend - nr_evap_tend) * dt);
  }

  qv.set(context, qv + qr2qv_evap_tend *dt);

  th_atm.set(context, th_atm + inv_exner*(-qr2qv_evap_tend * latvap * INV_CP) * dt);
}

} // namespace p3
} // namespace scream

#endif
