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
  //obtain minimum TKE allowed
  static constexpr auto mintke   = SC::mintke; // units:m2/s2

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

      //take max(mintke,tke(k))
      tke(k).set(tke(k) < mintke, mintke);

  });
}

} // namespace shoc
} // namespace scream

#endif
