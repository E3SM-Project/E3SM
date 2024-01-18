#ifndef P3_GET_LATENT_HEAT_IMPL_HPP
#define P3_GET_LATENT_HEAT_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
void Functions<S,D>
::get_latent_heat(const Int& nj, const Int& nk, view_2d<Spack>& v, view_2d<Spack>& s, view_2d<Spack>& f)
{
  constexpr Scalar latvap  = C::LatVap;
  constexpr Scalar latice  = C::LatIce;

  Kokkos::deep_copy(v, latvap);
  Kokkos::deep_copy(s, latvap + latice);
  Kokkos::deep_copy(f, latice);
}

} // namespace p3
} // namespace scream

#endif // P3_GET_LATENT_HEAT_IMPL_HPP
