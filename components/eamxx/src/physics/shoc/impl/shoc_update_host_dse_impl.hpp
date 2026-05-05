#ifndef SHOC_UPDATE_HOST_DSE_IMPL_HPP
#define SHOC_UPDATE_HOST_DSE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_host_dse(
  const MemberType& team,
  const Int& nlev,
  const uview_1d<const Pack>& thlm,
  const uview_1d<const Pack>& shoc_ql,
  const uview_1d<const Pack>& inv_exner,
  const uview_1d<const Pack>& zt_grid,
  const Scalar& phis,
  const uview_1d<Pack>& host_dse)
{
  const Int nlev_pack = ekat::npack<Pack>(nlev);

  // Constants used
  const auto lcond = C::LatVap.value;
  const auto cp    = C::CP.value;
  const auto ggr   = C::gravit.value;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
      Pack temp = (thlm(k)/inv_exner(k))+(lcond/cp)*shoc_ql(k);
      host_dse(k) = cp*temp+ggr*zt_grid(k)+phis;
  });
}

} // namespace shoc
} // namespace scream

#endif
