#ifndef ZM_ENTROPY_IMPL_HPP
#define ZM_ENTROPY_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm entropy. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
Real Functions<S,D>::entropy(
  // Inputs
  const MemberType& team,
  const Real& tk,
  const Real& p,
  const Real& qtot)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
  return 0;
}

} // namespace zm
} // namespace scream

#endif
