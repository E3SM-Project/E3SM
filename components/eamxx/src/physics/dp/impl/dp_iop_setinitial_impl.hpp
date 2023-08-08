#ifndef DP_IOP_SETINITIAL_IMPL_HPP
#define DP_IOP_SETINITIAL_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_setinitial. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::iop_setinitial(const Int& nelemd, const uview_1d<element_t>& elem)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
