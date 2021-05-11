
#ifndef SHOC_PBLINTD_INIT_POT_IMPL_HPP
#define SHOC_PBLINTD_INIT_POT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_pblintd_init_pot(
    const MemberType& team, const Int& nlev,
    const view_1d<const Spack>& thl, const view_1d<const Spack>& ql, const view_1d<const Spack>& q,
    const view_1d<Spack>& thv)
{
   // Compute virtual potential temperature
   const auto lcond = C::LatVap;
   const auto cp    = C::Cpair;
   const auto eps   = C::ZVIR; 
   const auto one   = C::ONE; 

   const Int nlev_pack = ekat::npack<Spack>(nlev);

   Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const int& k) {
     const auto th = thl(k) + (lcond/cp)*ql(k);
     thv(k) = th * (one + eps*q(k) - ql(k));
   });
}

} // namespace shoc
} // namespace scream

#endif
