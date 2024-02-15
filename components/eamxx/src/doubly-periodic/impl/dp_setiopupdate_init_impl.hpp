#ifndef DP_SETIOPUPDATE_INIT_IMPL_HPP
#define DP_SETIOPUPDATE_INIT_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp setiopupdate_init. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::setiopupdate_init()
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
