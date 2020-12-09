#ifndef SHOC_PBLINTD_SURF_TEMP_IMPL_HPP
#define SHOC_PBLINTD_SURF_TEMP_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_surf_temp. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_surf_temp(const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<const Spack>& z, const uview_1d<const Spack>& ustar, const uview_1d<const Spack>& obklen, const uview_1d<const Spack>& kbfs, const uview_1d<const Spack>& thv, const uview_1d<Spack>& tlv, const uview_1d<Spack>& pblh, const uview_1d<bool>& check, const uview_1d<Spack>& rino)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
