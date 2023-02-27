#ifndef SHOC_GRID_IMPL_HPP
#define SHOC_GRID_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_grid. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

/*
 * Define the thickness arrays of each column, to be used for
 * finite differencing throughout the SHOC parameterization,
 * also define air density in SHOC
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_grid(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& pdel,
  const uview_1d<Spack>&       dz_zt,
  const uview_1d<Spack>&       dz_zi,
  const uview_1d<Spack>&       rho_zt)
{
  const auto ggr = C::gravit;

  const auto s_zi_grid = ekat::scalarize(zi_grid);
  const auto s_zt_grid = ekat::scalarize(zt_grid);
  const auto s_dz_zi   = ekat::scalarize(dz_zi);

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

    // Compute shifts
    Spack zi_grid_k, zi_grid_kp1, zt_grid_k, zt_grid_km1;

    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack_m1 = range_pack;
    // index for _km1 should never go below 0
    range_pack_m1.set(range_pack < 1, 1);
    // and index for kp1 should never go above nlevi-1
    range_pack.set(range_pack>=(nlevi-1),nlevi-2);

    ekat::index_and_shift< 1>(s_zi_grid, range_pack,    zi_grid_k, zi_grid_kp1);
    ekat::index_and_shift<-1>(s_zt_grid, range_pack_m1, zt_grid_k, zt_grid_km1);

    // Define thickness of the thermodynamic gridpoints
    dz_zt(k) = zi_grid_k - zi_grid_kp1;

    // Define thickness of the interface grid points
    dz_zi(k) = zt_grid_km1 - zt_grid_k;

    // Define the air density on the thermo grid
    rho_zt(k) = (1/ggr)*(pdel(k)/dz_zt(k));
  });

  team.team_barrier();
  // Set lower condition for dz_zi
  s_dz_zi(0) = 0;
  s_dz_zi(nlevi-1) = s_zt_grid(nlev-1);
}

} // namespace shoc
} // namespace scream

#endif
