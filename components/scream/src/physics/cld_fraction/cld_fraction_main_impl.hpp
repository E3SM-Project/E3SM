#ifndef CLD_FRACTION_MAIN_IMPL_HPP
#define CLD_FRACTION_MAIN_IMPL_HPP

#include "physics/cld_fraction/cld_fraction_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace cldfrac {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void CldFractionFunctions<S,D>
::main(
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

    calc_icefrac(team,nk,oqi,oaist);

    calc_totalfrac(team,nk,oalst,oaist,oast);
  });
  Kokkos::fence();
} // main
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void CldFractionFunctions<S,D>
::calc_icefrac(
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
} // calc_icefrac
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void CldFractionFunctions<S,D>
::calc_totalfrac(
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
      ast(k) = max(aist(k),alst(k));
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // calc_totalfrac
/*-----------------------------------------------------------------*/

} // namespace cldfrac
} // namespace scream

#endif // CLD_FRACTION_MAIN_IMPL_HPP
