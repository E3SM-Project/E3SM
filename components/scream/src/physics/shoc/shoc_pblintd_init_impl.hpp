#ifndef SHOC_PBLINTD_INIT_IMPL_HPP
#define SHOC_PBLINTD_INIT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_init. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_init(const Scalar& z, bool& check, Scalar& rino, Scalar& pblh)
{
  const auto zero = C::ZERO;

  check = true;
  rino  = zero;
  pblh  = z;
}

} // namespace shoc
} // namespace scream

#endif
