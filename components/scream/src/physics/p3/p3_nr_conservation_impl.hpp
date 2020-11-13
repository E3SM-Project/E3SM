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
void Functions<S,D>::nr_conservation(const Spack& nr, const Spack& ni2nr_melt_tend, const Spack& nr_ice_shed_tend, const Spack& ncshdc, const Spack& nc2nr_autoconv_tend, const Spack& dt, Spack& nr_collect_tend, Spack& nr2ni_immers_freeze_tend, Spack& nr_selfcollect_tend, Spack& nr_evap_tend)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace p3
} // namespace scream

#endif
