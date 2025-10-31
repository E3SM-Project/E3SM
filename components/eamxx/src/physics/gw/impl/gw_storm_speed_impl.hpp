#ifndef GW_GW_STORM_SPEED_IMPL_HPP
#define GW_GW_STORM_SPEED_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/util/eamxx_utils.hpp"

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_storm_speed. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_storm_speed(
  // Inputs
  const MemberType& team,
  const GwCommonInit& init,
  const GwConvectInit& cinit,
  const Int& pver,
  const Real& storm_speed_min,
  const uview_1d<const Real>& ubm,
  const Int& mini,
  const Int& maxi,
  // Outputs
  Int& storm_speed,
  Real& uh,
  Real& umin,
  Real& umax)
{
  storm_speed = Int(Kokkos::copysign(ekat::impl::max(std::abs(ubm(cinit.k_src_wind))-storm_speed_min, Real(0)), ubm(cinit.k_src_wind)));

  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, maxi, mini+1), [&] (const int k, Real& lsum) {
      lsum += ubm(k)/(mini-maxi+1);
    }, Kokkos::Sum<Real>(uh));

  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&] {
    uh -= storm_speed;

    // Limit uh to table range.
    uh = ekat::impl::min(uh, Real(cinit.maxuh));
    uh = ekat::impl::max(uh, Real(-cinit.maxuh));
  });

  // Speeds for critical level filtering.
  if (maxi > mini) {
    umin =  init.pgwv*init.dc;
    umax = -init.pgwv*init.dc;
  }
  else {
    using ResultType = Kokkos::MinMax<Real>::value_type;
    ResultType min_max_result;

    Kokkos::parallel_reduce(
      Kokkos::TeamVectorRange(team, maxi, mini+1), [&] (const int k, ResultType& update) {
        if (ubm(k) < update.min_val) update.min_val = ubm(k);
        if (ubm(k) > update.max_val) update.max_val = ubm(k);
      }, Kokkos::MinMax<Real>(min_max_result));

    team.team_barrier();

    umin = min_max_result.min_val;
    umax = min_max_result.max_val;
  }
}

} // namespace gw
} // namespace scream

#endif
