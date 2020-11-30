#ifndef SHOC_COMPUTE_SHR_PROD_IMPL_HPP
#define SHOC_COMPUTE_SHR_PROD_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc compute_shr_prod. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_shr_prod(const Int& nlevi, const Int& nlev, const Int& shcol, const uview_1d<const Spack>& dz_zi, const uview_1d<const Spack>& u_wind, const uview_1d<const Spack>& v_wind, const uview_1d<Spack>& sterm)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
