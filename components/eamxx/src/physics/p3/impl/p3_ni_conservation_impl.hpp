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
void Functions<S,D>::ni_conservation(const Pack& ni, const Pack& ni_nucleat_tend, const Pack& nr2ni_immers_freeze_tend, 
const Pack& nc2ni_immers_freeze_tend, const Pack& ncheti_cnt, const Pack& nicnt, const Pack& ninuc_cnt, const Real& dt,
 Pack& ni2nr_melt_tend, Pack& ni_sublim_tend, Pack& ni_selfcollect_tend, const bool& use_hetfrz_classnuc, const Mask& context)
{
  const auto sink_ni = (ni2nr_melt_tend + ni_sublim_tend + ni_selfcollect_tend)*dt;

  Pack source_ni;
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
