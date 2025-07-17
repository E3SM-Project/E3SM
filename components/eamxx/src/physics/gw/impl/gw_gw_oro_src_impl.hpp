#ifndef GW_GW_ORO_SRC_IMPL_HPP
#define GW_GW_ORO_SRC_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_oro_src. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_oro_src(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const uview_1d<const Spack>& u,
const uview_1d<const Spack>& v,
const uview_1d<const Spack>& t,
const uview_1d<const Spack>& sgh,
const uview_1d<const Spack>& pmid,
const uview_1d<const Spack>& pint,
const uview_1d<const Spack>& dpm,
const uview_1d<const Spack>& zm,
const uview_1d<const Spack>& nm,
// Outputs
const uview_1d<Int>& src_level,
const uview_1d<Int>& tend_level,
const uview_1d<Spack>& tau,
const uview_1d<Spack>& ubm,
const uview_1d<Spack>& ubi,
const uview_1d<Spack>& xv,
const uview_1d<Spack>& yv,
const uview_1d<Spack>& c)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
