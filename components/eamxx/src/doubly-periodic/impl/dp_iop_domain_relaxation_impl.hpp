#ifndef DP_IOP_DOMAIN_RELAXATION_IMPL_HPP
#define DP_IOP_DOMAIN_RELAXATION_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_domain_relaxation. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::iop_domain_relaxation(const Int& nelemd, const Int& np, const Int& nlev, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const uview_1d<Spack>& dp, const Int& nelemd_todo, const Int& np_todo, const Spack& dt)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
