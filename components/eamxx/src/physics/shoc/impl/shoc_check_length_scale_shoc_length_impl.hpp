#ifndef SHOC_CHECK_LENGTH_SCALE_SHOC_LENGTH_IMPL_HPP
#define SHOC_CHECK_LENGTH_SCALE_SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::check_length_scale_shoc_length(
  const MemberType&      team,
  const Int&             nlev,
  const Scalar&          dx,
  const Scalar&          dy,
  const uview_1d<Spack>& shoc_mix)
{
  const auto minlen = scream::shoc::Constants<Real>::minlen;

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Ensure shoc_mix is in the interval [minval, sqrt(host_dx*host_dy)]
    shoc_mix(k) = ekat::min(std::sqrt(dx*dy),
                            ekat::max(minlen, shoc_mix(k)));
  });
}

} // namespace shoc
} // namespace scream

#endif
