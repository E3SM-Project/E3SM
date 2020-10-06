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
  const Scalar&          host_dx,
  const Scalar&          host_dy,
  const uview_1d<Spack>& shoc_mix)
{
  const auto minlen = scream::shoc::Constants<Real>::minlen;
  const auto maxlen = scream::shoc::Constants<Real>::maxlen;

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
    // Ensure shoc_mix is in the interval [minval, maxval]
    shoc_mix(k).set(maxlen < shoc_mix(k), maxlen);
    shoc_mix(k).set(minlen > shoc_mix(k), minlen);

    const auto tmp_val = sqrt(host_dx*host_dy);
    shoc_mix(k).set(tmp_val < shoc_mix(k), tmp_val);
  });
}

} // namespace shoc
} // namespace scream

#endif
