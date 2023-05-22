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
  const uview_1d<const Spack>& thlm,
  const uview_1d<const Spack>& shoc_ql,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& zt_grid,
  const Scalar& phis,
  const uview_1d<Spack>& host_dse)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);

  // Constants used
  const auto lcond = C::LatVap;
  const auto cp = C::CP;
  const auto ggr = C::gravit;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
      Spack temp = (thlm(k)/inv_exner(k))+(lcond/cp)*shoc_ql(k);
      host_dse(k) = cp*temp+ggr*zt_grid(k)+phis;
  });
}

} // namespace shoc
} // namespace scream

#endif
