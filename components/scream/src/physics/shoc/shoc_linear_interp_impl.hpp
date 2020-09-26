#ifndef SHOC_LINEAR_INTERP_IMPL_HPP
#define SHOC_LINEAR_INTERP_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc linear_interp. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::linear_interp(const uview_1d<const Spack>& x1, const uview_1d<const Spack>& x2, const uview_1d<const Spack>& y1, const uview_1d<Spack>& y2, const Int& km1, const Int& km2, const Int& ncol, const Spack& minthresh)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
