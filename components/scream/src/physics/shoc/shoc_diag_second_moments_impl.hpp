#ifndef SHOC_DIAG_SECOND_MOMENTS_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_second_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_second_moments(const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind, const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy, const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi, const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix, const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec, const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
