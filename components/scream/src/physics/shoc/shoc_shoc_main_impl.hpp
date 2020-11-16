#ifndef SHOC_SHOC_MAIN_IMPL_HPP
#define SHOC_SHOC_MAIN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_main. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_main(const Int& shcol, const Int& nlev, const Int& nlevi, const Spack& dtime, const Int& nadv, const uview_1d<const Spack>& host_dx, const uview_1d<const Spack>& host_dy, const uview_1d<const Spack>& thv, const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& pres, const uview_1d<const Spack>& presi, const uview_1d<const Spack>& pdel, const uview_1d<const Spack>& wthl_sfc, const uview_1d<const Spack>& wqw_sfc, const uview_1d<const Spack>& uw_sfc, const uview_1d<const Spack>& vw_sfc, const uview_1d<const Spack>& wtracer_sfc, const Int& num_qtracers, const uview_1d<const Spack>& w_field, const uview_1d<const Spack>& exner, const uview_1d<const Spack>& phis, const uview_1d<Spack>& host_dse, const uview_1d<Spack>& tke, const uview_1d<Spack>& thetal, const uview_1d<Spack>& qw, const uview_1d<Spack>& u_wind, const uview_1d<Spack>& v_wind, const uview_1d<Spack>& qtracers, const uview_1d<Spack>& wthv_sec, const uview_1d<Spack>& tkh, const uview_1d<Spack>& tk, const uview_1d<Spack>& shoc_ql, const uview_1d<Spack>& shoc_cldfrac, const uview_1d<Spack>& pblh, const uview_1d<Spack>& shoc_mix, const uview_1d<Spack>& isotropy, const uview_1d<Spack>& w_sec, const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& qwthl_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& w3, const uview_1d<Spack>& wqls_sec, const uview_1d<Spack>& brunt, const uview_1d<Spack>& shoc_ql2)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
