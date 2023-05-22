#ifndef P3_FIND_IMPL_HPP
#define P3_FIND_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 find functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
Int Functions<S,D>
::find_bottom (
    const MemberType& team,
    const uview_1d<const Scalar>& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present)
{
  log_present = false;
  Int k_xbot = 0;
  if (team.team_size() == 1) {
    for (Int k = kbot; k != ktop + kdir; k += kdir) {
      if (v(k) < small) continue;
      k_xbot = k;
      log_present = true;
      break;
    }
  } else {
    if (kdir == -1) {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, kbot - ktop + 1), [&] (Int k_, int& lmax) {
          const Int k = ktop + k_;
          if (v(k) >= small && k > lmax)
            lmax = k;
        }, Kokkos::Max<int>(k_xbot));
      log_present = k_xbot >= ktop;
    } else {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, ktop - kbot + 1), [&] (Int k_, int& lmin) {
          const Int k = kbot + k_;
          if (v(k) >= small && k < lmin)
            lmin = k;
        }, Kokkos::Min<int>(k_xbot));
      log_present = k_xbot <= ktop;
    }
  }
  return k_xbot;
}

template <typename S, typename D>
KOKKOS_FUNCTION
Int Functions<S,D>
::find_top (
    const MemberType& team,
    const uview_1d<const Scalar>& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present)
{
  log_present = false;
  Int k_xtop = 0;
  if (team.team_size() == 1) {
    for (Int k = ktop; k != kbot - kdir; k -= kdir) {
      if (v(k) < small) continue;
      k_xtop = k;
      log_present = true;
      break;
    }
  } else {
    if (kdir == -1) {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, kbot - ktop + 1), [&] (Int k_, int& lmin) {
          const Int k = ktop + k_;
          if (v(k) >= small && k < lmin)
            lmin = k;
        }, Kokkos::Min<int>(k_xtop));
      log_present = k_xtop <= kbot;
    } else {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, ktop - kbot + 1), [&] (Int k_, int& lmax) {
          const Int k = kbot + k_;
          if (v(k) >= small && k > lmax)
            lmax = k;
        }, Kokkos::Max<int>(k_xtop));
      log_present = k_xtop >= kbot;
    }
  }
  return k_xtop;
}

} // namespace p3
} // namespace scream

#endif
