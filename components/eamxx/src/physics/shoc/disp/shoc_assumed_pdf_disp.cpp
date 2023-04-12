#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_assumed_pdf_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const Int&                  nlevi,
  const view_2d<const Spack>& thetal,
  const view_2d<const Spack>& qw,
  const view_2d<const Spack>& w_field,
  const view_2d<const Spack>& thl_sec,
  const view_2d<const Spack>& qw_sec,
  const view_2d<const Spack>& wthl_sec,
  const view_2d<const Spack>& w_sec,
  const view_2d<const Spack>& wqw_sec,
  const view_2d<const Spack>& qwthl_sec,
  const view_2d<const Spack>& w3,
  const view_2d<const Spack>& pres,
  const view_2d<const Spack>& zt_grid,
  const view_2d<const Spack>& zi_grid,
  const WorkspaceMgr&         workspace_mgr,
  const view_2d<Spack>&       shoc_cldfrac,
  const view_2d<Spack>&       shoc_ql,
  const view_2d<Spack>&       wqls,
  const view_2d<Spack>&       wthv_sec,
  const view_2d<Spack>&       shoc_ql2)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    shoc_assumed_pdf(
      team, nlev, nlevi,
      ekat::subview(thetal, i),
      ekat::subview(qw, i),
      ekat::subview(w_field, i),
      ekat::subview(thl_sec, i),
      ekat::subview(qw_sec, i),
      ekat::subview(wthl_sec, i),
      ekat::subview(w_sec, i),
      ekat::subview(wqw_sec, i),
      ekat::subview(qwthl_sec, i),
      ekat::subview(w3, i),
      ekat::subview(pres, i),
      ekat::subview(zt_grid, i),
      ekat::subview(zi_grid, i),
      workspace,
      ekat::subview(shoc_cldfrac, i),
      ekat::subview(shoc_ql, i),
      ekat::subview(wqls, i),
      ekat::subview(wthv_sec, i),
      ekat::subview(shoc_ql2, i));
  });
}

} // namespace shoc
} // namespace scream
