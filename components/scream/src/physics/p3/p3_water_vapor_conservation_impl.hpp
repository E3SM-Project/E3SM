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
void Functions<S,D>::water_vapor_conservation(const Spack& qv, Spack& qidep, Spack& qinuc, const Spack& qi2qv_sublim_tend, const Spack& qr2qv_evap_tend, const Real& dt)
{
  const auto qv_avail = qv + (qi2qv_sublim_tend+qr2qv_evap_tend)*dt;
  const auto qv_sink  = (qidep + qinuc)*dt;

  const auto active = qv_sink > qv_avail && qv_sink > 1.e-20;
  if (active.any()) {
    const auto ratio = qv_avail / qv_sink;
    qidep.set(active, qidep*ratio);
    qinuc.set(active, qinuc*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
