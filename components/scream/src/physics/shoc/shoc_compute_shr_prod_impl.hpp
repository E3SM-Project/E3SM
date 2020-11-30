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
  const MemberType&            team  ,
  const Int&                   nlevi ,
  const Int&                   nlev  ,
  const Int&                   shcol ,
  const uview_1d<const Spack>& dz_zi ,
  const uview_1d<const Spack>& u_wind,
  const uview_1d<const Spack>& v_wind,
  const uview_1d<Spack>&       sterm)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);

  const auto sclr_uwind = scalarize(u_wind);//
  const auto sclr_vwind = scalarize(v_wind);

  //compute shear production term
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {

    auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack2 = range_pack1;
    range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway

    const Spack grid_dz = 1/dz_zi(k);
    // calculate vertical gradient of u&v wind
    Spack u_up_grid, u_grid, v_up_grid, v_grid;
    ekat::index_and_shift<-1>(sclr_uwind, range_pack2, u_grid, u_up_grid);
    ekat::index_and_shift<-1>(sclr_vwind, range_pack2, v_grid, v_up_grid);
    //const Spcak u_grad = grid_dz*(u_wind(km1)-u_wind(k));
    //const Spack v_grad = grid_dz*(v_wind(km1)-v_wind(k));
    //const Spack u_grad = grid_dz*(u_wind(k)-u_wind(k)+1);
    //const Spack v_grad = grid_dz*(v_wind(k)-v_wind(k)+1);
    const Spack u_grad(range_pack1 > 0 && range_pack1 < nlev, grid_dz*(u_up_grid - u_grid));
    const Spack v_grad(range_pack1 > 0 && range_pack1 < nlev, grid_dz*(v_up_grid - v_grid));
    sterm(k)           = u_grad*u_grad+v_grad*v_grad;
  });
  /*
   * Set lower and upper boundary for shear production
   * Note that the lower bound for shear production has already
   * been taken into account for the TKE boundary condition,
   * thus zero out here
    */
  sterm(0)[0]     = 0;
  const Int nlevi_pack = ekat::npack<Spack>(nlevi);
  const Int last_pack_entry = (nlevi%Spack::n == 0 ? Spack::n-1 : nlevi%Spack::n-1);
  sterm(nlevi_pack-1)[last_pack_entry] = 0;

}

} // namespace shoc
} // namespace scream

#endif
