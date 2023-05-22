#ifndef SHOC_COMPUTE_SHOC_VAPOR_IMPL_HPP
#define SHOC_COMPUTE_SHOC_VAPOR_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc compute_shoc_vapor. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * This function computes water vapor
 * based on SHOC's prognostic total water mixing ratio
 * and diagnostic cloud water mixing ratio.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_shoc_vapor(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& qw,
  const uview_1d<const Spack>& ql,
  const uview_1d<Spack>&       qv)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    qv(k) = qw(k) - ql(k);
  });
}

} // namespace shoc
} // namespace scream

#endif
