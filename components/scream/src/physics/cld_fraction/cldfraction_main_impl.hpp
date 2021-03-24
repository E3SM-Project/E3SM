#ifndef CLDFRACTION_MAIN_IMPL_HPP
#define CLDFRACTION_MAIN_IMPL_HPP

#include "physics/cld_fraction/cldfraction_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace cldfrac {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void Functions<S,D>
::cldfraction_main(
  const Int nj,
  const Int nk,
  const view_2d<const Spack>& qi,
  const view_2d<const Spack>& alst,
  const view_2d<Spack>& aist,
  const view_2d<Spack>& ast)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  Kokkos::parallel_for(
    "cld fraction main loop",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();

    const auto oqi   = ekat::subview(qi,   i);
    const auto oalst = ekat::subview(alst, i);
    const auto oaist = ekat::subview(aist, i);
    const auto oast  = ekat::subview(ast,  i);

    cldfraction_calc_icefrac(team,nk,oqi,oaist);

    cldfraction_calc_totalfrac(team,nk,oalst,oaist,oast);
  });
  Kokkos::fence();
} // cldfraction_main
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cldfraction_calc_icefrac(
  const MemberType& team,
  const Int& nk,
  const uview_1d<const Spack>& qi,
  const uview_1d<Spack>&       aist)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
      auto icecld = qi(k) > 1e-5;
      aist(k) = 0.0;
      aist(k).set(icecld, 1.0);
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // cldfraction_calc_icefrac
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cldfraction_calc_totalfrac(
  const MemberType& team,
  const Int& nk,
  const uview_1d<const Spack>& alst,
  const uview_1d<const Spack>& aist,
  const uview_1d<Spack>&       ast)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
      auto icecld = aist(k) > alst(k);
      ast(k).set(Smask(true), alst(k));
      ast(k).set(icecld, aist(k));
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // cldfraction_calc_totalfrac
/*-----------------------------------------------------------------*/

} // namespace cldfrac
} // namespace scream

#endif // CLDFRACTION_MAIN_IMPL_HPP
