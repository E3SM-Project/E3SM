#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_energy_integrals_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const view_2d<const Spack>& host_dse,
  const view_2d<const Spack>& pdel,
  const view_2d<const Spack>& rtm,
  const view_2d<const Spack>& rcm,
  const uview_2d<const Spack>& u_wind,
  const uview_2d<const Spack>& v_wind,
  const view_1d<Scalar>& se_b,
  const view_1d<Scalar>& ke_b,
  const view_1d<Scalar>& wv_b,
  const view_1d<Scalar>& wl_b)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    shoc_energy_integrals(team, nlev,
                          ekat::subview(host_dse, i),
                          ekat::subview(pdel, i),
                          ekat::subview(rtm, i),
                          ekat::subview(rcm, i),
                          ekat::subview(u_wind, i),
                          ekat::subview(v_wind, i),
                          se_b(i),
                          ke_b(i),
                          wv_b(i),
                          wl_b(i));
  });
}

} // namespace shoc
} // namespace scream
