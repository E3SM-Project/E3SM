#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_tke_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const view_2d<const Spack>&  wthv_sec,
  const view_2d<const Spack>&  shoc_mix,
  const view_2d<const Spack>&  dz_zi,
  const view_2d<const Spack>&  dz_zt,
  const view_2d<const Spack>&  pres,
  const view_2d<const Spack>&  u_wind,
  const view_2d<const Spack>&  v_wind,
  const view_2d<const Spack>&  brunt,
  const view_1d<const Scalar>& obklen,
  const view_2d<const Spack>&  zt_grid,
  const view_2d<const Spack>&  zi_grid,
  const view_1d<const Scalar>& pblh,
  const WorkspaceMgr&          workspace_mgr,
  const view_2d<Spack>&        tke,
  const view_2d<Spack>&        tk,
  const view_2d<Spack>&        tkh,
  const view_2d<Spack>&        isotropy)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);

    shoc_tke(team, nlev, nlevi, dtime,
             ekat::subview(wthv_sec, i),
             ekat::subview(shoc_mix, i),
             ekat::subview(dz_zi, i),
             ekat::subview(dz_zt, i),
             ekat::subview(pres, i),
             ekat::subview(u_wind, i),
             ekat::subview(v_wind, i),
             ekat::subview(brunt, i),
             obklen(i),
             ekat::subview(zt_grid, i),
             ekat::subview(zi_grid, i),
             pblh(i),
             workspace,
             ekat::subview(tke, i),
             ekat::subview(tk, i),
             ekat::subview(tkh, i),
             ekat::subview(isotropy, i));
  });
}

} // namespace shoc
} // namespace scream
