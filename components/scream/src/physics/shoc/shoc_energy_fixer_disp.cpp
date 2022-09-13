#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_energy_fixer_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const Int&                   nadv,
  const view_2d<const Spack>&  zt_grid,
  const view_2d<const Spack>&  zi_grid,
  const Int&                   se_b_slot,
  const Int&                   ke_b_slot,
  const Int&                   wv_b_slot,
  const Int&                   wl_b_slot,
  const Int&                   se_a_slot,
  const Int&                   ke_a_slot,
  const Int&                   wv_a_slot,
  const Int&                   wl_a_slot,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const Int&                   rho_zt_slot,
  const view_2d<const Spack>&  tke,
  const view_2d<const Spack>&  pint,
  const WorkspaceMgr&          workspace_mgr,
  const ScalarWorkspaceMgr&    workspace_mgr_local,
  const view_2d<Spack>&        host_dse)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);
    auto workspace_local = workspace_mgr_local.get_workspace(team);

    uview_1d<Scalar> se_b_int, ke_b_int, wv_b_int, wl_b_int, se_a_int, ke_a_int, wv_a_int, wl_a_int;
    se_b_int = workspace_local.get_space_in_slot(se_b_slot);
    ke_b_int = workspace_local.get_space_in_slot(ke_b_slot);
    wv_b_int = workspace_local.get_space_in_slot(wv_b_slot);
    wl_b_int = workspace_local.get_space_in_slot(wl_b_slot);
    se_a_int = workspace_local.get_space_in_slot(se_a_slot);
    ke_a_int = workspace_local.get_space_in_slot(ke_a_slot);
    wv_a_int = workspace_local.get_space_in_slot(wv_a_slot);
    wl_a_int = workspace_local.get_space_in_slot(wl_a_slot);

    uview_1d<Spack> rho_zt;
    rho_zt   = workspace.get_space_in_slot(rho_zt_slot);

    shoc_energy_fixer(team, nlev, nlevi, dtime, nadv,
                      ekat::subview(zt_grid, i),
                      ekat::subview(zi_grid, i),
                      se_b_int(i),
                      ke_b_int(i),
                      wv_b_int(i),
                      wl_b_int(i),
                      se_a_int(i),
                      ke_a_int(i),
                      wv_a_int(i),
                      wl_a_int(i),
                      wthl_sfc(i),
                      wqw_sfc(i),
                      rho_zt,
                      ekat::subview(tke, i),
                      ekat::subview(pint, i),
                      workspace,
                      ekat::subview(host_dse, i));

  });
}

} // namespace shoc
} // namespace scream
