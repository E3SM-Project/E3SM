#ifndef SHOC_COMPUTE_VERTICAL_SHEAR_TERMS_IMPL_HPP
#define SHOC_COMPUTE_VERTICAL_SHEAR_TERMS_IMPL_HPP

#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_vertical_shear_terms(
  const MemberType&              team,
  const Int&                     nlev,
  const Int&                     nlevi,
  const uview_1d<const Pack>&    dz_zi,
  const uview_1d<const Pack>&    u_wind,
  const uview_1d<const Pack>&    v_wind,
  const uview_1d<const Pack>&    w_field,
  const uview_1d<const Pack>&    zt_grid,
  const uview_1d<const Pack>&    zi_grid,
  const Workspace&               workspace,
  const uview_1d<Pack>&          du_dz_m,
  const uview_1d<Pack>&          dv_dz_m,
  const uview_1d<Pack>&          dw_dz_m)
{
  // Compute the SHOC-column vertical gradients on interfaces, then
  // interpolate them back to midpoint levels.
  uview_1d<Pack> du_dz_i, dv_dz_i, dw_dz_i;
  workspace.template take_many_contiguous_unsafe<3>(
    {"du_dz_i", "dv_dz_i", "dw_dz_i"},
    {&du_dz_i, &dv_dz_i, &dw_dz_i});

  const Int nlev_pack = ekat::npack<Pack>(nlev);
  const Int nlevi_pack = ekat::npack<Pack>(nlevi);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevi_pack), [&] (const Int& k) {
    du_dz_i(k) = 0;
    dv_dz_i(k) = 0;
    dw_dz_i(k) = 0;
  });
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    du_dz_m(k) = 0;
    dv_dz_m(k) = 0;
    dw_dz_m(k) = 0;
  });
  team.team_barrier();

  const auto s_u_wind  = scalarize(u_wind);
  const auto s_v_wind  = scalarize(v_wind);
  const auto s_w_field = scalarize(w_field);

  // Form the vertical gradients on the interface grid first so they are
  // consistent with SHOC's native staggered-grid treatment of shear production.
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    auto range_pack = ekat::range<IntPack>(k*Pack::n);
    const auto active_range = range_pack > 0 && range_pack < nlev;

    if (active_range.any()) {
      const Pack inv_dz = 1 / dz_zi(k);

      auto range_pack_safe = range_pack;
      range_pack_safe.set(range_pack < 1, 1);

      Pack u_grid, u_up_grid;
      Pack v_grid, v_up_grid;
      Pack w_grid, w_up_grid;
      ekat::index_and_shift<-1>(s_u_wind,  range_pack_safe, u_grid, u_up_grid);
      ekat::index_and_shift<-1>(s_v_wind,  range_pack_safe, v_grid, v_up_grid);
      ekat::index_and_shift<-1>(s_w_field, range_pack_safe, w_grid, w_up_grid);

      du_dz_i(k).set(active_range, inv_dz*(u_up_grid - u_grid));
      dv_dz_i(k).set(active_range, inv_dz*(v_up_grid - v_grid));
      dw_dz_i(k).set(active_range, inv_dz*(w_up_grid - w_grid));
    }
  });

  auto s_du_dz_i = scalarize(du_dz_i);
  auto s_dv_dz_i = scalarize(dv_dz_i);
  auto s_dw_dz_i = scalarize(dw_dz_i);
  s_du_dz_i(0) = 0;
  s_dv_dz_i(0) = 0;
  s_dw_dz_i(0) = 0;
  s_du_dz_i(nlevi-1) = 0;
  s_dv_dz_i(nlevi-1) = 0;
  s_dw_dz_i(nlevi-1) = 0;

  // Interpolate the interface-grid vertical gradients back to midpoint levels,
  // where SHOC carries thermodynamic variables and TKE.
  team.team_barrier();
  linear_interp(team, zi_grid, zt_grid, du_dz_i, du_dz_m, nlevi, nlev, 0);
  linear_interp(team, zi_grid, zt_grid, dv_dz_i, dv_dz_m, nlevi, nlev, 0);
  linear_interp(team, zi_grid, zt_grid, dw_dz_i, dw_dz_m, nlevi, nlev, 0);

  workspace.template release_many_contiguous<3>(
    {&du_dz_i, &dv_dz_i, &dw_dz_i});
}

} // namespace shoc
} // namespace scream

#endif
