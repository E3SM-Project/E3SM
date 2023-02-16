#ifndef SHOC_COMPUTE_BRUNT_SHOC_LENGTH_IMPL_HPP
#define SHOC_COMPUTE_BRUNT_SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_brunt_shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& thv,
  const uview_1d<const Spack>& thv_zi,
  const uview_1d<Spack>&       brunt)
{
  const auto ggr = C::gravit;
  const auto s_thv_zi = scalarize(thv_zi);

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Calculate thv_zi shift
    Spack thv_zi_k, thv_zi_kp1;
    auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack2 = range_pack1;
    range_pack2.set(range_pack1 >= (nlevi-1), nlevi-2);
    ekat::index_and_shift<1>(s_thv_zi, range_pack2, thv_zi_k, thv_zi_kp1);

    brunt(k) = (ggr/thv(k))*(thv_zi_k-thv_zi_kp1)/dz_zt(k);
  });
}

} // namespace shoc
} // namespace scream

#endif
