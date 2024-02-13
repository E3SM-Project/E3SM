#ifndef DP_APPLY_IOP_FORCING_IMPL_HPP
#define DP_APPLY_IOP_FORCING_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp apply_iop_forcing. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::apply_iop_forcing(const Int& nelemd, const uview_1d<element_t>& elem, hvcoord_t& hvcoord, const hybrid_t& hybrid, const timelevel_t& tl, const Int& n, const bool& t_before_advance, const Int& nets, const Int& nete)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
