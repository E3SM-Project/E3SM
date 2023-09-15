#ifndef DP_CRM_RESOLVED_TURB_IMPL_HPP
#define DP_CRM_RESOLVED_TURB_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp crm_resolved_turb. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::crm_resolved_turb(const Int& nelemd, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const Int& nelemd_todo, const Int& np_todo)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
