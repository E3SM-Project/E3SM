#ifndef P3_NI_CONSERVATION_IMPL_HPP
#define P3_NI_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ni_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::ni_conservation(const Spack& ni, const Spack& ni_nucleat_tend, const Spack& nr2ni_immers_freeze_tend, 
const Spack& nc2ni_immers_freeze_tend, const Spack& ncheti_cnt, const Spack& nicnt, const Spack& ninuc_cnt, const Real& dt,
 Spack& ni2nr_melt_tend, Spack& ni_sublim_tend, Spack& ni_selfcollect_tend, const bool& use_hetfrz_classnuc, const Smask& context)
{
  const auto sink_ni = (ni2nr_melt_tend + ni_sublim_tend + ni_selfcollect_tend)*dt;

  Spack source_ni;
  if(use_hetfrz_classnuc){
    source_ni = ni + (ni_nucleat_tend+nr2ni_immers_freeze_tend+ncheti_cnt+nicnt+ninuc_cnt)*dt;
  }
   else {
    source_ni = ni + (ni_nucleat_tend+nr2ni_immers_freeze_tend+nc2ni_immers_freeze_tend)*dt;
   }
  const auto mask = sink_ni > source_ni && context;
  if (mask.any()) {
    const auto ratio = source_ni/sink_ni;
    ni2nr_melt_tend.set(mask, ni2nr_melt_tend*ratio);
    ni_sublim_tend.set(mask, ni_sublim_tend*ratio);
    ni_selfcollect_tend.set(mask, ni_selfcollect_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
