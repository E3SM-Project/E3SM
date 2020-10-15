#ifndef SHOC_COMPUTE_SHOC_VAPOR_IMPL_HPP
#define SHOC_COMPUTE_SHOC_VAPOR_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc compute_shoc_vapor. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_shoc_vapor(const Int& shcol, const Int& nlev, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& ql, const uview_1d<Spack>& qv)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
