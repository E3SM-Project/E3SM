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
  const view_1d<const Scalar>& se_b,
  const view_1d<const Scalar>& ke_b,
  const view_1d<const Scalar>& wv_b,
  const view_1d<const Scalar>& wl_b,
  const view_1d<const Scalar>& se_a,
  const view_1d<const Scalar>& ke_a,
  const view_1d<const Scalar>& wv_a,
  const view_1d<const Scalar>& wl_a,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const view_2d<const Spack>&  rho_zt,
  const view_2d<const Spack>&  tke,
  const view_2d<const Spack>&  pint,
  const WorkspaceMgr&          workspace_mgr,
  const view_2d<Spack>&        host_dse)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);

    shoc_energy_fixer(team, nlev, nlevi, dtime, nadv,
                      ekat::subview(zt_grid, i),
                      ekat::subview(zi_grid, i),
                      se_b(i),
                      ke_b(i),
                      wv_b(i),
                      wl_b(i),
                      se_a(i),
                      ke_a(i),
                      wv_a(i),
                      wl_a(i),
                      wthl_sfc(i),
                      wqw_sfc(i),
                      ekat::subview(rho_zt, i),
                      ekat::subview(tke, i),
                      ekat::subview(pint, i),
                      workspace,
                      ekat::subview(host_dse, i));

  });
}

} // namespace shoc
} // namespace scream
