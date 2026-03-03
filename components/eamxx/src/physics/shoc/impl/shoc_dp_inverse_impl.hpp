#ifndef SHOC_DP_INVERSE_IMPL_HPP
#define SHOC_DP_INVERSE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::dp_inverse(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Pack>& rho_zt,
  const uview_1d<const Pack>& dz_zt,
  const uview_1d<Pack>&       rdp_zt)
{
  const auto ggr = C::gravit.value;

  const Int nlev_pack = ekat::npack<Pack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    rdp_zt(k) = 1/(ggr*rho_zt(k)*dz_zt(k));
  });
}

} // namespace shoc
} // namespace scream

#endif
