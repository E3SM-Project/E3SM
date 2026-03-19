
#ifndef SHOC_PBLINTD_INIT_POT_IMPL_HPP
#define SHOC_PBLINTD_INIT_POT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "share/physics/physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_pblintd_init_pot(
    const MemberType& team, const Int& nlev,
    const view_1d<const Pack>& thl, const view_1d<const Pack>& ql, const view_1d<const Pack>& q,
    const view_1d<Pack>& thv)
{
   // Compute virtual potential temperature
   const auto lcond = C::LatVap.value;
   const auto cp    = C::Cpair.value;
   const auto eps   = C::ZVIR;
   const auto one   = C::ONE;

   const Int nlev_pack = ekat::npack<Pack>(nlev);

   Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const int& k) {
     const auto th = thl(k) + (lcond/cp)*ql(k);
     thv(k) = th * (one + eps*q(k) - ql(k));
   });
}

} // namespace shoc
} // namespace scream

#endif
