#include "shoc_functions.hpp"

#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::diag_third_shoc_moments_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const Int&                  nlevi,
  const Scalar&               c_diag_3rd_mom,
  const bool&                 shoc_1p5tke,
  const view_2d<const Pack>& w_sec,
  const view_2d<const Pack>& thl_sec,
  const view_2d<const Pack>& wthl_sec,
  const view_2d<const Pack>& isotropy,
  const view_2d<const Pack>& brunt,
  const view_2d<const Pack>& thetal,
  const view_2d<const Pack>& tke,
  const view_2d<const Pack>& dz_zt,
  const view_2d<const Pack>& dz_zi,
  const view_2d<const Pack>& zt_grid,
  const view_2d<const Pack>& zi_grid,
  const WorkspaceMgr&         workspace_mgr,
  const view_2d<Pack>&       w3)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const auto nlev_packs = ekat::npack<Pack>(nlev);
  const auto policy = TPF::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    diag_third_shoc_moments(
      team, nlev, nlevi,
      c_diag_3rd_mom,
      shoc_1p5tke,
      ekat::subview(w_sec, i),
      ekat::subview(thl_sec, i),
      ekat::subview(wthl_sec, i),
      ekat::subview(isotropy, i),
      ekat::subview(brunt, i),
      ekat::subview(thetal, i),
      ekat::subview(tke, i),
      ekat::subview(dz_zt, i),
      ekat::subview(dz_zi, i),
      ekat::subview(zt_grid, i),
      ekat::subview(zi_grid, i),
      workspace,
      ekat::subview(w3, i));
  });
}

} // namespace shoc
} // namespace scream
