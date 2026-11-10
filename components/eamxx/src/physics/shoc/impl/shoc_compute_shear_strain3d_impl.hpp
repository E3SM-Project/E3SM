#ifndef SHOC_COMPUTE_SHEAR_STRAIN3D_IMPL_HPP
#define SHOC_COMPUTE_SHEAR_STRAIN3D_IMPL_HPP

#include "shoc_assemble_shoc_shear_strain3d_impl.hpp"
#include "shoc_compute_vertical_shear_terms_impl.hpp"

namespace scream {
namespace shoc {

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
    uview_1d<Pack> du_dz_m, dv_dz_m, dw_dz_m;
    workspace.template take_many_contiguous_unsafe<3>(
      {"du_dz_m", "dv_dz_m", "dw_dz_m"},
      {&du_dz_m, &dv_dz_m, &dw_dz_m});

    compute_vertical_shear_terms(team, nlev, nlevi,
                                 Kokkos::subview(dz_zi, i, Kokkos::ALL()),
                                 Kokkos::subview(u_wind, i, Kokkos::ALL()),
                                 Kokkos::subview(v_wind, i, Kokkos::ALL()),
                                 Kokkos::subview(w_field, i, Kokkos::ALL()),
                                 Kokkos::subview(zt_grid, i, Kokkos::ALL()),
                                 Kokkos::subview(zi_grid, i, Kokkos::ALL()),
                                 workspace,
                                 du_dz_m, dv_dz_m, dw_dz_m);

    assemble_shoc_shear_strain3d(team, nlev,
                                 Kokkos::subview(shear_strain3d_components, i, Kokkos::ALL(), Kokkos::ALL()),
                                 du_dz_m, dv_dz_m, dw_dz_m,
                                 Kokkos::subview(shear_strain3d, i, Kokkos::ALL()));

    workspace.template release_many_contiguous<3>(
      {&du_dz_m, &dv_dz_m, &dw_dz_m});
  });
}
#endif

} // namespace shoc
} // namespace scream

#endif
