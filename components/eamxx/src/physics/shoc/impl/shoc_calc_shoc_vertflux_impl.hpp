#ifndef SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP
#define SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_shoc_vertflux(
  const MemberType& team,
  const Int& nlev,
  const uview_1d<const Pack>& tkh_zi,
  const uview_1d<const Pack>& dz_zi,
  const uview_1d<const Pack>& invar,
  const uview_1d<Pack>& vertflux)
{
  const Int nlev_pack = ekat::npack<Pack>(nlev);
  const auto sinvar = scalarize(invar);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    auto range_pack1 = ekat::range<IntSmallPack>(k*Pack::n);
    auto range_pack2 = range_pack1;
    range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway

    Pack up_grid, grid;
    ekat::index_and_shift<-1>(sinvar, range_pack2, grid, up_grid);
    const auto active_range = range_pack1 > 0 && range_pack1 < nlev;
    if (active_range.any()) {
      const Pack grid_dz = 1 / dz_zi(k); // vertical grid diff squared
      // Compute the vertical flux via downgradient diffusion
      vertflux(k).set(active_range, -(tkh_zi(k) * grid_dz * (up_grid - grid)));
    }
  });
}

} // namespace shoc
} // namespace scream

#endif
