#ifndef SHOC_CALC_SHOC_VARORCOVAR_IMPL_HPP
#define SHOC_CALC_SHOC_VARORCOVAR_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_shoc_varorcovar(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar& tunefac,
  const uview_1d<const Spack>& isotropy_zi,
  const uview_1d<const Spack>& tkh_zi,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& invar1,
  const uview_1d<const Spack>& invar2,
  const uview_1d<Spack>&       varorcovar)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto sinvar1 = scalarize(invar1);
  const auto sinvar2 = scalarize(invar2);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Calculate shift
    Spack invar1_s, invar1_sm1, invar2_s, invar2_sm1;
    auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack2 = range_pack1;
    range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway
    ekat::index_and_shift<-1>(sinvar1, range_pack2, invar1_s, invar1_sm1);
    ekat::index_and_shift<-1>(sinvar2, range_pack2, invar2_s, invar2_sm1);

    const auto active_range = range_pack1 > 0 && range_pack1 < nlev;
    if (active_range.any()) {
      const Spack grid_dz = 1 / dz_zi(k);
      const Spack grid_dz2 = ekat::square(grid_dz); // vertical grid diff squared

      // Compute the variance or covariance
      varorcovar(k).set(active_range, tunefac*(isotropy_zi(k)*tkh_zi(k))*grid_dz2*(invar1_sm1 - invar1_s)*(invar2_sm1 - invar2_s));
    }
  });
}

} // namespace shoc
} // namespace scream

#endif
