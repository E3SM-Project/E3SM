#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::pblintd_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   npbl,
  const view_2d<const Spack>&  z,
  const view_2d<const Spack>&  zi,
  const view_2d<const Spack>&  thl,
  const view_2d<const Spack>&  ql,
  const view_2d<const Spack>&  q,
  const view_2d<const Spack>&  u,
  const view_2d<const Spack>&  v,
  const view_1d<const Scalar>& ustar,
  const view_1d<const Scalar>& obklen,
  const view_1d<const Scalar>& kbfs,
  const view_2d<const Spack>&  cldn,
  const WorkspaceMgr&          workspace_mgr,
  const view_1d<Scalar>&       pblh)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace       = workspace_mgr.get_workspace(team);

    Scalar pblh_;
    pblintd(team, nlev, nlevi, npbl,
            ekat::subview(z, i),
            ekat::subview(zi, i),
            ekat::subview(thl, i),
            ekat::subview(ql, i),
            ekat::subview(q, i),
            ekat::subview(u, i),
            ekat::subview(v, i),
            ustar(i),
            obklen(i),
            kbfs(i),
            ekat::subview(cldn, i),
            workspace,
            pblh_);
    pblh(i) = pblh_;
  });
}

} // namespace shoc
} // namespace scream
