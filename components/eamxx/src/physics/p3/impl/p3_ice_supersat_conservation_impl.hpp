#ifndef P3_ICE_SUPERSAT_CONSERVATION_IMPL_HPP
#define P3_ICE_SUPERSAT_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice_supersat_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::ice_supersat_conservation(Spack& qv2qi_vapdep_tend, Spack& qv2qi_nucleat_tend, Spack& qinuc_cnt, const Spack& cld_frac_i, const Spack& qv, const Spack& qv_sat_i, const Spack& t_atm, const Real& dt, const Spack& qi2qv_sublim_tend, const Spack& qr2qv_evap_tend, const bool& use_hetfrz_classnuc, const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar cp     = C::CP;
  constexpr Scalar rv     = C::RH2O;
  constexpr Scalar latvap = C::LatVap;
  constexpr Scalar latice = C::LatIce;
  constexpr Scalar latsublim2 = (latvap+latice)*(latvap+latice);

  Spack qv_sink;
  if(use_hetfrz_classnuc){
      qv_sink = qv2qi_vapdep_tend + qv2qi_nucleat_tend + qinuc_cnt; // in [kg/kg] cell-avg values
  }
  else{
      qv_sink = qv2qi_vapdep_tend + qv2qi_nucleat_tend; // in [kg/kg] cell-avg values
  }

  const auto mask = qv_sink > qsmall && cld_frac_i > 1e-20 && context;
  if (mask.any()) {
    // --- Available water vapor for deposition/nucleation
    auto qv_avail = (qv + (qi2qv_sublim_tend+qr2qv_evap_tend)*dt - qv_sat_i) /
      (1 + latsublim2*qv_sat_i / (cp*rv*square(t_atm)) ) / dt;

    // --- Only excess water vapor can be limited
    qv_avail = max(qv_avail, 0);

    const auto sink_gt_avail = qv_sink > qv_avail && mask;
    if (sink_gt_avail.any()) {
      const auto fract = qv_avail / qv_sink;
      if(use_hetfrz_classnuc){
         qinuc_cnt.set(sink_gt_avail, qinuc_cnt * fract);
      }
      qv2qi_nucleat_tend.set(sink_gt_avail, qv2qi_nucleat_tend * fract);
      qv2qi_vapdep_tend.set(sink_gt_avail, qv2qi_vapdep_tend * fract);
    }
  }
}

} // namespace p3
} // namespace scream

#endif
