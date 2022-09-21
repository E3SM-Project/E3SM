#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::compute_shoc_vapor_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const view_2d<const Spack>& qw,
  const view_2d<const Spack>& ql,
  const view_2d<Spack>&       qv)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    compute_shoc_vapor(team, nlev,
                       ekat::subview(qw, i),
                       ekat::subview(ql, i),
                       ekat::subview(qv, i));
  });
}

} // namespace shoc
} // namespace scream
