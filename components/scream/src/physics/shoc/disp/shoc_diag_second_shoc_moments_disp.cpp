#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::diag_second_shoc_moments_disp(
  const Int& shcol, const Int& nlev, const Int& nlevi,
  const view_2d<const Spack>& thetal,
  const view_2d<const Spack>& qw,
  const view_2d<const Spack>& u_wind,
  const view_2d<const Spack>& v_wind,
  const view_2d<const Spack>& tke,
  const view_2d<const Spack>& isotropy,
  const view_2d<const Spack>& tkh,
  const view_2d<const Spack>& tk,
  const view_2d<const Spack>& dz_zi,
  const view_2d<const Spack>& zt_grid,
  const view_2d<const Spack>& zi_grid,
  const view_2d<const Spack>& shoc_mix,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const view_1d<const Scalar>& uw_sfc,
  const view_1d<const Scalar>& vw_sfc,
  const view_1d<Scalar>& ustar2,
  const view_1d<Scalar>& wstar,
  const WorkspaceMgr& workspace_mgr,
  const view_2d<Spack>& thl_sec,
  const view_2d<Spack>& qw_sec,
  const view_2d<Spack>& wthl_sec,
  const view_2d<Spack>& wqw_sec,
  const view_2d<Spack>& qwthl_sec,
  const view_2d<Spack>& uw_sec,
  const view_2d<Spack>& vw_sec,
  const view_2d<Spack>& wtke_sec,
  const view_2d<Spack>& w_sec)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    diag_second_shoc_moments(
      team, nlev, nlevi,
      ekat::subview(thetal, i),
      ekat::subview(qw, i),
      ekat::subview(u_wind, i),
      ekat::subview(v_wind, i),
      ekat::subview(tke, i),
      ekat::subview(isotropy, i),
      ekat::subview(tkh, i),
      ekat::subview(tk, i),
      ekat::subview(dz_zi, i),
      ekat::subview(zt_grid, i),
      ekat::subview(zi_grid, i),
      ekat::subview(shoc_mix, i),
      wthl_sfc(i),
      wqw_sfc(i),
      uw_sfc(i),
      vw_sfc(i),
      ustar2(i),
      wstar(i),
      workspace,
      ekat::subview(thl_sec, i),
      ekat::subview(qw_sec, i),
      ekat::subview(wthl_sec, i),
      ekat::subview(wqw_sec, i),
      ekat::subview(qwthl_sec, i),
      ekat::subview(uw_sec, i),
      ekat::subview(vw_sec, i),
      ekat::subview(wtke_sec, i),
      ekat::subview(w_sec, i));
  });
}

} // namespace shoc
} // namespace scream
