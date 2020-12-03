#ifndef SHOC_PBLINTD_IMPL_HPP
#define SHOC_PBLINTD_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd(const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<const Spack>& z, const uview_1d<const Spack>& zi, const uview_1d<const Spack>& thl, const uview_1d<const Spack>& ql, const uview_1d<const Spack>& q, const uview_1d<const Spack>& u, const uview_1d<const Spack>& v, const uview_1d<const Spack>& ustar, const uview_1d<const Spack>& obklen, const uview_1d<const Spack>& kbfs, const uview_1d<const Spack>& cldn, const uview_1d<Spack>& pblh)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
