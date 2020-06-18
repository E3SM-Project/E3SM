#ifndef P3_FUNCTIONS_UPDATE_PROGNOSTICS_IMPL_HPP
#define P3_FUNCTIONS_UPDATE_PROGNOSTICS_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_prognostic_ice(
  const Spack& qcheti,      const Spack& qccol,  const Spack& qcshd,     const Spack& nccol,
  const Spack& ncheti,      const Spack& ncshdc, const Spack& qrcol,     const Spack& nrcol,
  const Spack& qrheti,      const Spack& nrheti, const Spack& nrshdr,    const Spack& qimlt,
  const Spack& nimlt,       const Spack& qisub,  const Spack& qidep,     const Spack& qinuc,
  const Spack& ninuc,       const Spack& nislf,  const Spack& nisub,     const Spack& qiberg,
  const Spack& exner,       const Spack& xxls,   const Spack& xlf,       const bool log_predictNc,
  const Smask& log_wetgrowth, const Scalar dt,  const Scalar& nmltratio, const Spack& rhorime_c,
  Spack& th, Spack& qv, Spack& qitot, Spack& nitot, Spack& qirim, Spack& birim, Spack& qc,
  Spack& nc, Spack& qr, Spack& nr,
  const Smask& context)
{
  constexpr Scalar QSMALL          = C::QSMALL;
  constexpr Scalar INV_RHO_RIMEMAX = C::INV_RHO_RIMEMAX;

  qc.set(context, qc + (-qcheti-qccol-qcshd-qiberg)*dt);
  if ( log_predictNc ){
    nc.set(context, nc + (-nccol-ncheti)*dt);
  }

  qr.set(context, qr + (-qrcol+qimlt-qrheti+qcshd)*dt);

  //apply factor to source for rain number from melting of ice, (ad-hoc
  // but accounts for rapid evaporation of small melting ice particles)
  nr.set(context, nr + (-nrcol-nrheti+nmltratio*nimlt+nrshdr+ncshdc)*dt);

  const auto qitot_not_small = qitot >= QSMALL && context;

  if ( qitot_not_small.any() ) {
    birim.set(qitot_not_small, birim - ((qisub + qimlt) / qitot) * dt * birim);
    qirim.set(qitot_not_small, qirim - ((qisub + qimlt) * qirim / qitot) * dt);
    qitot.set(qitot_not_small, qitot - (qisub + qimlt) * dt);
  }

  const auto dum = (qrcol + qccol + qrheti + qcheti) * dt;
  qitot.set(context, qitot + (qidep + qinuc + qiberg) * dt + dum);
  qirim.set(context, qirim + dum);

  birim.set(context, birim + (qrcol * INV_RHO_RIMEMAX + qccol / rhorime_c + (qrheti +
                                                                             qcheti) * INV_RHO_RIMEMAX) * dt);

  nitot.set(context, nitot + (ninuc - nimlt - nisub - nislf + nrheti + ncheti) * dt);

  //PMC nCat deleted interactions_loop

  const auto qirim_lt_thresh = qirim < 0 && context;
  if (qirim_lt_thresh.any()){
    qirim.set(qirim_lt_thresh, 0);
    birim.set(qirim_lt_thresh, 0);
  }

  // densify under wet growth
  // -- to be removed post-v2.1.  Densification automatically happens
  //    during wet growth due to parameterized rime density --

  qirim.set(log_wetgrowth && context, qitot);
  birim.set(log_wetgrowth && context, qirim * INV_RHO_RIMEMAX);

  // densify in above freezing conditions and melting
  // -- future work --
  //   Ideally, this will be treated with the predicted liquid fraction in ice.
  //   Alternatively, it can be simplified by tending qirim -- qitot
  //   and birim such that rho_rim (qirim/birim) --> rho_liq during melting.
  // ==
  qv.set(context, qv + (-qidep+qisub-qinuc)*dt);

  constexpr Scalar INV_CP = C::INV_CP;
  th.set(context, th + exner * ((qidep - qisub + qinuc) * xxls * INV_CP +
                                (qrcol + qccol + qcheti + qrheti - qimlt + qiberg) * xlf * INV_CP) * dt);
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_prognostic_liquid(
  const Spack& qcacc, const Spack& ncacc,
  const Spack& qcaut,const Spack& ncautc, const Spack& ncautr,
  const Spack& ncslf, const Spack& qrevp, const Spack& nrevp, const Spack& nrslf,
  const bool log_predictNc, const Spack& inv_rho, const Spack& exner, const Spack& xxlv,
  const Scalar dt, Spack& th, Spack& qv, Spack& qc, Spack& nc, Spack& qr, Spack& nr,
  const Smask& context)
{
  constexpr Scalar NCCNST = C::NCCNST;
  constexpr int IPARAM    = C::IPARAM;
  constexpr Scalar INV_CP = C::INV_CP;

  qc.set(context, qc + (-qcacc-qcaut)*dt);
  qr.set(context, qr + (qcacc+qcaut-qrevp)*dt);

  if (log_predictNc) {
    nc.set(context, nc + (-ncacc-ncautc+ncslf)*dt);
  }
  else {
    nc.set(context, NCCNST * inv_rho);
  }

  if (IPARAM == 1 || IPARAM == 2) {
    nr.set(context, nr + (sp(0.5) * ncautc - nrslf - nrevp) * dt);
  }
  else {
    nr.set(context, nr + (ncautr - nrslf - nrevp) * dt);
  }

  qv.set(context, qv + qrevp *dt);

  th.set(context, th + exner*(-qrevp * xxlv * INV_CP) * dt);
}

} // namespace p3
} // namespace scream

#endif
