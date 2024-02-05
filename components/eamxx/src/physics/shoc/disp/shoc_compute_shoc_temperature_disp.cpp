#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::compute_shoc_temperature_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const view_2d<const Spack>& thetal,
  const view_2d<const Spack>& ql,
  const view_2d<const Spack>& inv_exner,
  const view_2d<Spack>&       tabs)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    compute_shoc_temperature(team, nlev,
                             ekat::subview(thetal, i),
                             ekat::subview(ql, i),
		             ekat::subview(inv_exner, i),
                             ekat::subview(tabs, i));
  });
}

} // namespace shoc
} // namespace scream
