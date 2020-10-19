#ifndef SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP
#define SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc update_prognostics_implicit. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::update_prognostics_implicit(const Int& shcol, const Int& nlev, const Int& nlevi, const Int& num_tracer, const Spack& dtime, const uview_1d<const Spack>& dz_zt, const uview_1d<const Spack>& dz_zi, const uview_1d<const Spack>& rho_zt, const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& uw_sfc, const uview_1d<const Spack>& vw_sfc, const uview_1d<const Spack>& wthl_sfc, const uview_1d<const Spack>& wqw_sfc, const uview_1d<Spack>& thetal, const uview_1d<Spack>& qw, const uview_1d<Spack>& tracer, const uview_1d<Spack>& tke, const uview_1d<Spack>& u_wind, const uview_1d<Spack>& v_wind)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
