#ifndef SHOC_CLIPPING_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP
#define SHOC_CLIPPING_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::clipping_diag_third_shoc_moments(
  const MemberType& team,
  const Int& nlevi,
  const uview_1d<const Spack>& w_sec_zi,
  const uview_1d<Spack>& w3)
{
  const Int nlevi_pack = ekat::npack<Spack>(nlevi);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevi_pack), [&] (const Int& k) {
    const auto w3clip = scream::shoc::Constants<Scalar>::w3clip;
    Spack tsign(1);

    const auto theterm = w_sec_zi(k);
    const auto cond    = w3clip*ekat::sqrt(2*ekat::cube(theterm));
    const auto w3clipdef = 0.02;

    tsign.set(w3(k)<0, -1);
    w3(k).set(tsign*w3(k) > cond, w3clipdef);
  });
}

} // namespace shoc
} // namespace scream

#endif
