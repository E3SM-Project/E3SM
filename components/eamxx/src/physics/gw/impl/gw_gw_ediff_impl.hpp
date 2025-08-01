#ifndef GW_GW_EDIFF_IMPL_HPP
#define GW_GW_EDIFF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_ediff. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_ediff(
// Inputs
const Int& ncol,
const Int& pver,
const Int& kbot,
const Int& ktop,
const uview_1d<const Int>& tend_level,
const uview_1d<const Spack>& gwut,
const uview_1d<const Spack>& ubm,
const uview_1d<const Spack>& nm,
const uview_1d<const Spack>& rho,
const Spack& dt,
const Spack& gravit,
const uview_1d<const Spack>& pmid,
const uview_1d<const Spack>& rdpm,
const uview_1d<const Spack>& c,
// Outputs
const uview_1d<Spack>& egwdffi)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
