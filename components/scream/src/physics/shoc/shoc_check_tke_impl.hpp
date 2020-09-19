#ifndef SHOC_CHECK_TKE_IMPL_HPP
#define SHOC_CHECK_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::check_tke(
  const MemberType& team,
  const Int& nlev,
  const uview_1d<Spack>& tke)
{
  const Int nlev_pack = ekat::pack::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {

      //tke(k).set(tke(k)<0.1, mintke);

  });
}

} // namespace shoc
} // namespace scream

#endif
