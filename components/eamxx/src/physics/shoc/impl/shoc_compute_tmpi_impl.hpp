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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    const auto mask  = ekat::range<IntSmallPack>(k*Spack::n) > 0;
    if (mask.any()) {
      tmpi(k).set(mask, dtime*(ggr*rho_zi(k))/dz_zi(k));
    }
  });
}

} // namespace shoc
} // namespace scream

#endif
