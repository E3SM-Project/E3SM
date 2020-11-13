#ifndef P3_WATER_VAPOR_CONSERVATION_IMPL_HPP
#define P3_WATER_VAPOR_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 water_vapor_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::water_vapor_conservation(const Spack& qv, Spack& qidep, Spack& qinuc, const Spack& qi2qv_sublim_tend, const Spack& qr2qv_evap_tend, const Spack& dt)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace p3
} // namespace scream

#endif
