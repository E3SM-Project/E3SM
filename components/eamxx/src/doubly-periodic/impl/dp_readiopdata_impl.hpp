#ifndef DP_READIOPDATA_IMPL_HPP
#define DP_READIOPDATA_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp readiopdata. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::readiopdata(const Int& plev, const bool& iop_update_phase1, const uview_1d<const Spack>& hyam, const uview_1d<const Spack>& hybm)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
