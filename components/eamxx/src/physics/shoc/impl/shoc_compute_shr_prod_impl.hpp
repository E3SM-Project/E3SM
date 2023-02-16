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
void Functions<S,D>
::compute_shr_prod(
  const MemberType&            team,
  const Int&                   nlevi,
  const Int&                   nlev,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& u_wind,
  const uview_1d<const Spack>& v_wind,
  const uview_1d<Spack>&       sterm)
{
  // Turbulent coefficient
  const Scalar Ck_sh = 0.1;

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  //scalarize so that we can use shift to compute the differece  ( x(k-1) - x )
  const auto sclr_uwind = scalarize(u_wind); //for u_wind
  const auto sclr_vwind = scalarize(v_wind); //for v_wind

  //compute shear production term
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

    auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
    const auto active_range = range_pack1 > 0 && range_pack1 < nlev;
    if (active_range.any()) {
      const Spack grid_dz = 1/dz_zi(k);

      // calculate vertical gradient of u&v wind
      auto range_pack2 = range_pack1;
      range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway
      Spack u_up_grid, u_grid, v_up_grid, v_grid;
      ekat::index_and_shift<-1>(sclr_uwind, range_pack2, u_grid, u_up_grid); //for u_wind
      ekat::index_and_shift<-1>(sclr_vwind, range_pack2, v_grid, v_up_grid); //for v_wind
      const Spack u_grad(active_range, grid_dz*(u_up_grid - u_grid));
      const Spack v_grad(active_range, grid_dz*(v_up_grid - v_grid));

      //compute shear production
      sterm(k).set(active_range, Ck_sh*(u_grad*u_grad+v_grad*v_grad));
    }
  });
  /*
   * Set lower and upper boundary for shear production
   * Note that the lower bound for shear production has already
   * been taken into account for the TKE boundary condition,
   * thus zero out here
    */
  const auto s_sterm = ekat::scalarize(sterm);
  s_sterm(0)       = 0;
  s_sterm(nlevi-1) = 0;
}

} // namespace shoc
} // namespace scream

#endif
