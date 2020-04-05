#ifndef P3_FUNCTIONS_GET_LATENT_HEAT_IMPL_HPP
#define P3_FUNCTIONS_GET_LATENT_HEAT_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_functions_math_impl.hpp"


namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::get_latent_heat(const Int& ni, const Int& nk, const MemberType& team, view_2d<Spack>& v, view_2d<Spack>& s, view_2d<Spack>& f) 
{
   constexpr Scalar lapvap  = C::LatVap;
   constexpr Scalar latice  = C::LatIce;
  
   Kokkos::parallel_for (Kokkos::TeamThreadRange (team, ni), [=] (int& i) {
     Kokkos::parallel_for (Kokkos::ThreadVectorRange(team, nk), [=] (int& k) {
        v(i,k) = lapvap;
        s(i,k) = lapvap + latice;
        f(i,k) = latice;
     });
   });
}

} // namespace p3
} // namespace scream

#endif
