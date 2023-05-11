#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_length_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const view_1d<const Scalar>& dx,
  const view_1d<const Scalar>& dy,
  const view_2d<const Spack>&  zt_grid,
  const view_2d<const Spack>&  zi_grid,
  const view_2d<const Spack>&  dz_zt,
  const view_2d<const Spack>&  tke,
  const view_2d<const Spack>&  thv,
  const WorkspaceMgr&          workspace_mgr,
  const view_2d<Spack>&        brunt,
  const view_2d<Spack>&        shoc_mix)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);

    shoc_length(team, nlev, nlevi, dx(i), dy(i),
                ekat::subview(zt_grid, i),
                ekat::subview(zi_grid, i),
                ekat::subview(dz_zt, i),
                ekat::subview(tke, i),
                ekat::subview(thv, i),
                workspace,
                ekat::subview(brunt, i),
                ekat::subview(shoc_mix, i));
  });
}

} // namespace shoc
} // namespace scream
