#ifndef ZM_IENTROPY_IMPL_HPP
#define ZM_IENTROPY_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm ientropy. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::ientropy(
  // Inputs
  const MemberType& team,
  const Int& rcall,
  const Real& s,
  const Real& p,
  const Real& qt,
  const Real& tfg,
  // Outputs
  Real& t,
  Real& qst)
{
  // TODO
}

} // namespace zm
} // namespace scream

#endif
