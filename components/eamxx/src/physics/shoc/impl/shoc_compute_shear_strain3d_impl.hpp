#ifndef SHOC_COMPUTE_SHEAR_STRAIN3D_IMPL_HPP
#define SHOC_COMPUTE_SHEAR_STRAIN3D_IMPL_HPP

#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_shear_strain3d(
  const MemberType&              team,
  const Int&                     nlev,
  const Int&                     nlevi,
  const uview_2d<const Pack>&    shear_strain3d_components,
  const uview_1d<const Pack>&    dz_zi,
  const uview_1d<const Pack>&    u_wind,
  const uview_1d<const Pack>&    v_wind,
  const uview_1d<const Pack>&    w_field,
  const uview_1d<const Pack>&    zt_grid,
  const uview_1d<const Pack>&    zi_grid,
  const Workspace&               workspace,
  const uview_1d<Pack>&          shear_strain3d)
{
  uview_1d<Pack> du_dz_i, dv_dz_i, dw_dz_i, du_dz_m, dv_dz_m, dw_dz_m;
  workspace.template take_many_contiguous_unsafe<6>(
    {"du_dz_i", "dv_dz_i", "dw_dz_i", "du_dz_m", "dv_dz_m", "dw_dz_m"},
    {&du_dz_i, &dv_dz_i, &dw_dz_i, &du_dz_m, &dv_dz_m, &dw_dz_m});

  const Int nlev_pack = ekat::npack<Pack>(nlev);
  const Int nlevi_pack = ekat::npack<Pack>(nlevi);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevi_pack), [&] (const Int& k) {
    du_dz_i(k) = 0;
    dv_dz_i(k) = 0;
    dw_dz_i(k) = 0;
    du_dz_m(k) = 0;
    dv_dz_m(k) = 0;
    dw_dz_m(k) = 0;
  });
  team.team_barrier();

  const auto s_u_wind  = scalarize(u_wind);
  const auto s_v_wind  = scalarize(v_wind);
  const auto s_w_field = scalarize(w_field);

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

  team.team_barrier();
  linear_interp(team, zi_grid, zt_grid, du_dz_i, du_dz_m, nlevi, nlev, 0);
  linear_interp(team, zi_grid, zt_grid, dv_dz_i, dv_dz_m, nlevi, nlev, 0);
  linear_interp(team, zi_grid, zt_grid, dw_dz_i, dw_dz_m, nlevi, nlev, 0);

  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    const Pack A00 = shear_strain3d_components(0,k);
    const Pack A01 = shear_strain3d_components(1,k);
    const Pack A10 = shear_strain3d_components(2,k);
    const Pack A11 = shear_strain3d_components(3,k);
    const Pack A20 = shear_strain3d_components(4,k);
    const Pack A21 = shear_strain3d_components(5,k);

    const Pack A02 = du_dz_m(k);
    const Pack A12 = dv_dz_m(k);
    const Pack A22 = dw_dz_m(k);

    const Pack S00 = A00;
    const Pack S11 = A11;
    const Pack S22 = A22;
    const Pack S01 = 0.5 * (A01 + A10);
    const Pack S02 = 0.5 * (A02 + A20);
    const Pack S12 = 0.5 * (A12 + A21);

    shear_strain3d(k) =
      2.0 * (S00*S00 + S11*S11 + S22*S22
           + 2.0*S01*S01 + 2.0*S02*S02 + 2.0*S12*S12);
  });

  workspace.template release_many_contiguous<6>(
    {&du_dz_i, &dv_dz_i, &dw_dz_i, &du_dz_m, &dv_dz_m, &dw_dz_m});
}

#ifdef SCREAM_SHOC_SMALL_KERNELS
template<typename S, typename D>
void Functions<S,D>::compute_shear_strain3d_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const Int&                  nlevi,
  const view_3d<const Pack>&  shear_strain3d_components,
  const view_2d<const Pack>&  dz_zi,
  const view_2d<const Pack>&  u_wind,
  const view_2d<const Pack>&  v_wind,
  const view_2d<const Pack>&  w_field,
  const view_2d<const Pack>&  zt_grid,
  const view_2d<const Pack>&  zi_grid,
  const WorkspaceMgr&         workspace_mgr,
  const view_2d<Pack>&        shear_strain3d)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const auto nlev_packs = ekat::npack<Pack>(nlev);
  const auto policy = TPF::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    auto workspace = workspace_mgr.get_workspace(team);

    compute_shear_strain3d(team, nlev, nlevi,
                           Kokkos::subview(shear_strain3d_components, i, Kokkos::ALL(), Kokkos::ALL()),
                           Kokkos::subview(dz_zi, i, Kokkos::ALL()),
                           Kokkos::subview(u_wind, i, Kokkos::ALL()),
                           Kokkos::subview(v_wind, i, Kokkos::ALL()),
                           Kokkos::subview(w_field, i, Kokkos::ALL()),
                           Kokkos::subview(zt_grid, i, Kokkos::ALL()),
                           Kokkos::subview(zi_grid, i, Kokkos::ALL()),
                           workspace,
                           Kokkos::subview(shear_strain3d, i, Kokkos::ALL()));
  });
}
#endif

} // namespace shoc
} // namespace scream

#endif
