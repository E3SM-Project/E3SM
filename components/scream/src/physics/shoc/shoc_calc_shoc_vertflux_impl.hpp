#ifndef SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP
#define SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
void Functions<S,D>
::calc_shoc_vertflux(
  const MemberType& team,
  const Int& shcol, const Int& nlev,
  const uview_2d<const Spack>& tkh_zi,
  const uview_2d<const Spack>& dz_zi,
  const uview_2d<const Spack>& invar,
  const uview_2d<Spack>& vertflux)
{
  const Int nlev_pack = scream::pack::npack<Spack>(nlev);
  const auto sinvar = scalarize(invar);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, shcol), [&] (const Int& i) {
    const auto osinvar = util::subview(sinvar, i);
    for (Int k = 0; k < nlev_pack; ++k) {
      auto range_pack1 = scream::pack::range<IntSmallPack>(k*Spack::n);
      auto range_pack2 = range_pack1;
      range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway

      Spack up_grid, grid;
      pack::index_and_shift<-1>(osinvar, range_pack2, grid, up_grid);
      const Spack grid_dz = 1 / dz_zi(i,k); // vertical grid diff squared
      // Compute the vertical flux via downgradient diffusion
      vertflux(i,k).set(range_pack1 > 0 && range_pack1 < nlev, -(tkh_zi(i,k) * grid_dz * (up_grid - grid)));
    };
  });
}

} // namespace shoc
} // namespace scream

#endif
