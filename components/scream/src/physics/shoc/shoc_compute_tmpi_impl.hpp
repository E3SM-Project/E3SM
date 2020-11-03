#ifndef SHOC_COMPUTE_TMPI_IMPL_HPP
#define SHOC_COMPUTE_TMPI_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_tmpi(
  const MemberType&            team,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const uview_1d<const Spack>& rho_zi,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<Spack>&       tmpi)
{
  const auto ggr = C::gravit;

  tmpi(0)[0] = 0;

  const Int nlev_pack = ekat::npack<Spack>(nlevi);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    tmpi(k).set(range_pack > 0, dtime*(ggr*rho_zi(k))/dz_zi(k));
  });
}

} // namespace shoc
} // namespace scream

#endif
