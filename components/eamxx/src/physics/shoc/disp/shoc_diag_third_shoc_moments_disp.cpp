#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::diag_third_shoc_moments_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const Int&                  nlevi,
  const view_2d<const Spack>& w_sec,
  const view_2d<const Spack>& thl_sec,
  const view_2d<const Spack>& wthl_sec,
  const view_2d<const Spack>& isotropy,
  const view_2d<const Spack>& brunt,
  const view_2d<const Spack>& thetal,
  const view_2d<const Spack>& tke,
  const view_2d<const Spack>& dz_zt,
  const view_2d<const Spack>& dz_zi,
  const view_2d<const Spack>& zt_grid,
  const view_2d<const Spack>& zi_grid,
  const WorkspaceMgr&         workspace_mgr,
  const view_2d<Spack>&       w3)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    diag_third_shoc_moments(
      team, nlev, nlevi,
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
