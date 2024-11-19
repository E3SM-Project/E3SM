

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

template <>
void Functions<Real, DefaultDevice>
::check_values_disp(const uview_2d<const Spack>& qv, const uview_2d<const Spack>& temp, const Int& ktop, const Int& kbot,
               const Int& timestepcount, const bool& force_abort, const Int& source_ind,
               const uview_2d<const Scalar>& col_loc, const Int& nj, const Int& nk)
{

  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);

  Kokkos::parallel_for(
    "p3_check_values",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();

    check_values(ekat::subview(qv, i), ekat::subview(temp, i), ktop, kbot, timestepcount, force_abort, 900,
                 team, ekat::subview(col_loc, i));

  });

}

} // namespace p3
} // namespace scream
