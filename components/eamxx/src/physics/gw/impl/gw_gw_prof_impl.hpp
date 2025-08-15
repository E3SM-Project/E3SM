#ifndef GW_GW_PROF_IMPL_HPP
#define GW_GW_PROF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_prof. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_prof(
// Inputs
const Int& pver,
const Int& ncol,
const Spack& cpair,
const uview_1d<const Spack>& t,
const uview_1d<const Spack>& pmid,
const uview_1d<const Spack>& pint,
// Outputs
const uview_1d<Spack>& rhoi,
const uview_1d<Spack>& ti,
const uview_1d<Spack>& nm,
const uview_1d<Spack>& ni)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
