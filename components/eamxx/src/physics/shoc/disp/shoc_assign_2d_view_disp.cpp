#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_assign_2d_view_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const view_2d<const Spack>&  input_view,
  const view_2d<Spack>&        output_view)
{
  /*using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const Int& k) {
        output_view (i,k) = input_view(i,k);
    });
  });*/
}

} // namespace shoc
} // namespace scream
