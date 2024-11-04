#ifndef SHOC_COMPUTE_SHOC_TEMPERATURE_IMPL_HPP
#define SHOC_COMPUTE_SHOC_TEMPERATURE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc compute_shoc_temperature. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * This function computes absolute temperature 
 * based on SHOC's prognostic liquid water potential temperature.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_shoc_temperature(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& thetal,
  const uview_1d<const Spack>& ql,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<Spack>&       tabs)
{

  const Scalar cp = C::CP;
  const Scalar lcond = C::LatVap;

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    tabs(k) = thetal(k)/inv_exner(k)+(lcond/cp)*ql(k);
  });
}

} // namespace shoc
} // namespace scream

#endif
