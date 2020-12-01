#ifndef P3_WATER_VAPOR_CONSERVATION_IMPL_HPP
#define P3_WATER_VAPOR_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 water_vapor_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::water_vapor_conservation(const Spack& qv, Spack& qv2qi_vapdep_tend, Spack& qv2qi_nucleat_tend, const Spack& qi2qv_sublim_tend, const Spack& qr2qv_evap_tend, const Real& dt, const Smask& context)
{
  const auto qv_avail = qv + (qi2qv_sublim_tend+qr2qv_evap_tend)*dt;
  const auto qv_sink  = (qv2qi_vapdep_tend + qv2qi_nucleat_tend)*dt;

  const auto mask = qv_sink > qv_avail && qv_sink > 1.e-20 && context;
  if (mask.any()) {
    const auto ratio = qv_avail / qv_sink;
    qv2qi_vapdep_tend.set(mask, qv2qi_vapdep_tend*ratio);
    qv2qi_nucleat_tend.set(mask, qv2qi_nucleat_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
