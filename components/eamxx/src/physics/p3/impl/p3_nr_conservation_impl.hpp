#ifndef P3_NR_CONSERVATION_IMPL_HPP
#define P3_NR_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 nr_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::nr_conservation(const Spack& nr, const Spack& ni2nr_melt_tend, const Spack& nr_ice_shed_tend, const Spack& ncshdc, const Spack& nc2nr_autoconv_tend, const Real& dt, const Real& nmltratio, Spack& nr_collect_tend, Spack& nr2ni_immers_freeze_tend, Spack& nr_selfcollect_tend, Spack& nr_evap_tend, const Smask& context)
{
  const auto sink_nr = (nr_collect_tend + nr2ni_immers_freeze_tend + nr_selfcollect_tend + nr_evap_tend)*dt;
  const auto source_nr = nr + (ni2nr_melt_tend*nmltratio + nr_ice_shed_tend + ncshdc + nc2nr_autoconv_tend)*dt;
  const auto mask = sink_nr > source_nr && context;
  if (mask.any()) {
    const auto ratio = source_nr/sink_nr;
    nr_collect_tend.set(mask, nr_collect_tend*ratio);
    nr2ni_immers_freeze_tend.set(mask, nr2ni_immers_freeze_tend*ratio);
    nr_selfcollect_tend.set(mask, nr_selfcollect_tend*ratio);
    nr_evap_tend.set(mask, nr_evap_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
