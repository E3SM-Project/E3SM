#include "shoc_functions.hpp"

#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_tke_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const Scalar&                lambda_low,
  const Scalar&                lambda_high,
  const Scalar&                lambda_slope,
  const Scalar&                lambda_thresh,
  const Scalar&                Ckh,
  const Scalar&                Ckm,
  const bool&                  shoc_1p5tke,
  const bool&                  do_3d_turb,
  const view_2d<const Pack>&  wthv_sec,
  const view_2d<const Pack>&  strain2,
  const view_2d<const Pack>&  shoc_mix,
  const view_2d<const Pack>&  dz_zi,
  const view_2d<const Pack>&  dz_zt,
  const view_2d<const Pack>&  pres,
  const view_2d<const Pack>&  tabs,
  const view_2d<const Pack>&  u_wind,
  const view_2d<const Pack>&  v_wind,
  const view_2d<const Pack>&  brunt,
  const view_2d<const Pack>&  zt_grid,
  const view_2d<const Pack>&  zi_grid,
  const view_1d<const Scalar>& pblh,
  const WorkspaceMgr&          workspace_mgr,
  const view_2d<Pack>&        tke,
  const view_2d<Pack>&        tk,
  const view_2d<Pack>&        tkh,
  const view_2d<Pack>&        isotropy)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const auto nlev_packs = ekat::npack<Pack>(nlev);
  const auto policy = TPF::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);

    shoc_tke(team, nlev, nlevi, dtime,
             lambda_low, lambda_high, lambda_slope, lambda_thresh,
             Ckh, Ckm, shoc_1p5tke, do_3d_turb,
             ekat::subview(wthv_sec, i),
             ekat::subview(strain2, i),
             ekat::subview(shoc_mix, i),
             ekat::subview(dz_zi, i),
             ekat::subview(dz_zt, i),
             ekat::subview(pres, i),
             ekat::subview(tabs, i),
             ekat::subview(u_wind, i),
             ekat::subview(v_wind, i),
             ekat::subview(brunt, i),
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
