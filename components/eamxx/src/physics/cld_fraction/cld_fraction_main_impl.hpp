#ifndef CLD_FRACTION_MAIN_IMPL_HPP
#define CLD_FRACTION_MAIN_IMPL_HPP

#include "physics/cld_fraction/cld_fraction_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace cld_fraction {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void CldFractionFunctions<S,D>
::main(
  const Int nj,
  const Int nk,
  const view_2d<const Spack>& qi,
  const view_2d<const Spack>& liq_cld_frac,
  const view_2d<Spack>& ice_cld_frac,
  const view_2d<Spack>& tot_cld_frac)
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
    const auto oliq_cld_frac = ekat::subview(liq_cld_frac, i);
    const auto oice_cld_frac = ekat::subview(ice_cld_frac, i);
    const auto otot_cld_frac = ekat::subview(tot_cld_frac,  i);

    calc_icefrac(team,nk,oqi,oice_cld_frac);

    calc_totalfrac(team,nk,oliq_cld_frac,oice_cld_frac,otot_cld_frac);
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
  const uview_1d<Spack>&       ice_cld_frac)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {
      const Real ice_frac_threshold = 1e-12;
      auto icecld = qi(k) > ice_frac_threshold;
      ice_cld_frac(k) = 0.0;
      ice_cld_frac(k).set(icecld, 1.0);
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
  const uview_1d<const Spack>& liq_cld_frac,
  const uview_1d<const Spack>& ice_cld_frac,
  const uview_1d<Spack>&       tot_cld_frac)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {
      tot_cld_frac(k) = max(ice_cld_frac(k),liq_cld_frac(k));
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // calc_totalfrac
/*-----------------------------------------------------------------*/

} // namespace cld_fraction
} // namespace scream

#endif // CLD_FRACTION_MAIN_IMPL_HPP
