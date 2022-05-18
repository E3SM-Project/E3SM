#ifndef P3_NC_CONSERVATION_IMPL_HPP
#define P3_NC_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 nc_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::nc_conservation(const Spack& nc, const Spack& nc_selfcollect_tend, const Real& dt, Spack& nc_collect_tend, Spack& nc2ni_immers_freeze_tend, Spack& nc_accret_tend, Spack& nc2nr_autoconv_tend, const Smask& context)
{
  const auto sink_nc = (nc_collect_tend + nc2ni_immers_freeze_tend + nc_accret_tend + nc2nr_autoconv_tend)*dt;
  const auto source_nc = nc + nc_selfcollect_tend*dt;
  const auto mask = sink_nc > source_nc && context;
  if (mask.any()) {
    const auto ratio = source_nc/sink_nc;
    nc_collect_tend.set(mask, nc_collect_tend*ratio);
    nc2ni_immers_freeze_tend.set(mask, nc2ni_immers_freeze_tend*ratio);
    nc_accret_tend.set(mask, nc_accret_tend*ratio);
    nc2nr_autoconv_tend.set(mask, nc2nr_autoconv_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
