#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::update_host_dse_disp(
  const Int& shcol,
  const Int& nlev,
  const view_2d<const Spack>& thlm,
  const view_2d<const Spack>& shoc_ql,
  const view_2d<const Spack>& inv_exner,
  const view_2d<const Spack>& zt_grid,
  const view_1d<const Scalar>& phis,
  const view_2d<Spack>& host_dse)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    update_host_dse(
      team, nlev,
      ekat::subview(thlm, i),
      ekat::subview(shoc_ql, i),
      ekat::subview(inv_exner, i),
      ekat::subview(zt_grid, i),
      phis(i),
      ekat::subview(host_dse, i));
  });
}

} // namespace shoc
} // namespace scream
