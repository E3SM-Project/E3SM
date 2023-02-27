#ifndef SHOC_LINEAR_INTERP_IMPL_HPP
#define SHOC_LINEAR_INTERP_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc linear_interp. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::linear_interp(
  const MemberType& team,
  const uview_1d<const Spack>& x1,
  const uview_1d<const Spack>& x2,
  const uview_1d<const Spack>& y1,
  const uview_1d<Spack>& y2,
  const Int& km1,
  const Int& km2,
  const Scalar& minthresh)
{
  const auto sx1 = scalarize(x1);
  const auto sy1 = scalarize(y1);
  const Int km2_pack = ekat::npack<Spack>(km2);

  if (km1 == km2+1) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, km2_pack), [&] (const Int& k2) {
      Spack x1, x1s, y1, y1s; // s->-1 shift
      auto indx_pack = ekat::range<IntSmallPack>(k2*Spack::n + 1);

      // Avoid reading out of bounds for x1/y1
      indx_pack.set(indx_pack>=km1, km1-1);

      ekat::index_and_shift<-1>(sx1, indx_pack, x1, x1s);
      ekat::index_and_shift<-1>(sy1, indx_pack, y1, y1s);

      y2(k2) = y1s + (y1-y1s)*(x2(k2)-x1s)/(x1-x1s);
    });
  }
  else if (km2 == km1+1) {

    // TODO: This implementation may need to be tweaked for better CPU performance.
    // Specifically, CPU perf may improve if the top and bottom pack are handled in
    // a Kokkos::single block so that all conditionals can be removed from the ||4.
    // Keep thread0 unoccupied so it can process the Kokkos::single while other threads
    // are doing the ||4.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, km2_pack), [&] (const Int& k2) {
      Spack x1, x1s, y1, y1s; // s->-1 shift
      auto indx_pack = ekat::range<IntSmallPack>(k2*Spack::n);
      indx_pack.set(indx_pack < 1, 1); // special shift for 0 boundary case
      indx_pack.set(indx_pack >= km2-1, km1-1); // special shift for top boundary case
      ekat::index_and_shift<-1>(sx1, indx_pack, x1, x1s);
      ekat::index_and_shift<-1>(sy1, indx_pack, y1, y1s);

      y2(k2) = y1s + (y1-y1s)*(x2(k2)-x1s)/(x1-x1s);
    });
  }
  else {
    EKAT_KERNEL_REQUIRE_MSG(false, "Unsupported dimensions for linear interp");
  }
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, km2_pack), [&] (const Int& k2) {
    y2(k2).set(y2(k2) < minthresh, minthresh);
  });
}

} // namespace shoc
} // namespace scream

#endif
