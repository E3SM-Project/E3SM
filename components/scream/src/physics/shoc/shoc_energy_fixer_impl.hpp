#ifndef SHOC_ENERGY_FIXER_IMPL_HPP
#define SHOC_ENERGY_FIXER_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_energy_fixer. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_energy_fixer(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const Int&                   nadv,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const Scalar&                se_b,
  const Scalar&                ke_b,
  const Scalar&                wv_b,
  const Scalar&                wl_b,
  const Scalar&                se_a,
  const Scalar&                ke_a,
  const Scalar&                wv_a,
  const Scalar&                wl_a,
  const Scalar&                wthl_sfc,
  const Scalar&                wqw_sfc,
  const uview_1d<const Spack>& rho_zt,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& pint,
  const uview_1d<Spack>&       rho_zi,
  const uview_1d<Spack>&       host_dse)
{
  // Constants
  const auto cp = C::CP;
  const auto lcond = C::LatVap;
  const auto lice = C::LatIce;
  const auto mintke = SC::mintke;
  const auto ggr = C::gravit;

  // Local variables
  Scalar te_a = 0;
  Scalar te_b = 0;
  Scalar se_dis = 0;
  Int shoctop = 0;

  // Compute linear interpolation of data into rho_zi
  linear_interp(team,zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,0);
  team.team_barrier();

  // Compute the host timestep
  const Scalar hdtime = dtime*float(nadv);

  // Compute the total energy before and after SHOC call
  const auto s_rho_zi = ekat::scalarize(rho_zi);
  const Scalar shf = wthl_sfc*cp*s_rho_zi(nlevi-1);
  const Scalar lhf = wqw_sfc*s_rho_zi(nlevi-1);
  te_a = se_a + ke_a + (lcond+lice)*wv_a +lice*wl_a;
  te_b = se_b + ke_b + (lcond+lice)*wv_b + lice*wl_b;
  te_b += (shf+lhf*(lcond+lice))*hdtime;

  // Limit the energy fixer to find highest layer where SHOC is active.
  // Find first level where tke is higher than lowest threshold.
  const auto s_tke = ekat::scalarize(tke);
  const auto s_pint = ekat::scalarize(pint);
  while (s_tke(shoctop) == mintke && shoctop < nlev-2) {
    shoctop += 1;
  }

  // Compute the disbalance of total energy, over depth where SHOC is active.
  se_dis = (te_a - te_b)/(s_pint(nlevi-1) - s_pint(shoctop));

  // Update host_dse
  const int begin_pack = shoctop/Spack::n;
  const auto nlev_packs = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, begin_pack, nlev_packs), [&] (const Int& k) {
    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);

    host_dse(k).set(range_pack >= shoctop && range_pack < nlev, host_dse(k)-se_dis*ggr);
  });
}

} // namespace shoc
} // namespace scream

#endif
