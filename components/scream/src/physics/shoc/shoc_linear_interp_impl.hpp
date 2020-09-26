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
  const Int& ncol,
  const Scalar& minthresh)
{
  const auto sx1 = scalarize(x1);
  const auto sx2 = scalarize(x2);
  const auto sy1 = scalarize(y1);
  const auto sy2 = scalarize(y2);

  if (km1 == km2+1) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, km2), [&] (const Int& k2) {
      const Int k1 = k2+1;
      sy2(k2) = sy1(k1-1) + (sy1(k1)-sy1(k1-1))*(sx2(k2)-sx1(k1-1))/(sx1(k1)-sx1(k1-1));
    });
  }
  else if (km2 == km1+1) {
    sy2(0) = sy1(0) + (sy1(1)-sy1(0))*(sx2(0)-sx1(0))/(sx1(1)-sx1(0));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, km2-2), [&] (const Int& k2) {
      const Int k2s = k2+1;
      const Int k1 = k2s;
      sy2(k2s) = sy1(k1-1) + (sy1(k1)-sy1(k1-1))*(sx2(k2s)-sx1(k1-1))/(sx1(k1)-sx1(k1-1));
    });
    const Int k2 = km2-1;
    const Int k1 = km1-1;
    sy2(k2) = sy1(k1) + (sy1(k1)-sy1(k1-1))*(sx2(k2)-sx1(k1))/(sx1(k1)-sx1(k1-1));
  }
  else {
    EKAT_KERNEL_REQUIRE_MSG(false, "Unsupported dimensions for linear interp");
  }
  team.team_barrier();

  const Int km2_pack = ekat::npack<Spack>(km2);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, km2_pack), [&] (const Int& k2) {
    y2(k2).set(y2(k2) < minthresh, minthresh);
  });
}

} // namespace shoc
} // namespace scream

#endif
