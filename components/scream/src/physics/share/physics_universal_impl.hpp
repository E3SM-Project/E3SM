#ifndef PHYSICS_UNIVERSAL_IMPL_HPP
#define PHYSICS_UNIVERSAL_IMPL_HPP

#include "physics_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace physics {

using PC = Constants<Real>;
/*
 * Implementation of universal physics functions. Clients should NOT #include
 * this file, #include physics_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_exner(const Spack& pmid, const Smask& range_mask)
{
  Spack result;
  const Spack exner = 1.0 / pow( pmid/PC::P0, PC::RD*PC::INV_CP );
  result.set(range_mask,exner);
  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::T_to_th(const Spack& T_atm, const Spack& exner, const Smask& range_mask)
{
  Spack result;
  const Spack th_atm = T_atm*exner;
  result.set(range_mask,th_atm);
  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::th_to_T(const Spack& th_atm, const Spack& exner, const Smask& range_mask)
{
  Spack result;
  const Spack T_atm = th_atm/exner;
  result.set(range_mask,T_atm);
  return result;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_dz(const Spack& zi_top, const Spack& zi_bot, const Smask& range_mask)
{
  Spack result;
  const Spack dz = zi_top-zi_bot;
  result.set(range_mask,dz);
  return result;
}


} // namespace physics
} // namespace scream

#endif // PHYSICS_UNIVERSAL_IMPL_HPP
