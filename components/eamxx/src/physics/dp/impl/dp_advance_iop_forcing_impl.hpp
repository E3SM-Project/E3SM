#ifndef DP_ADVANCE_IOP_FORCING_IMPL_HPP
#define DP_ADVANCE_IOP_FORCING_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp advance_iop_forcing. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::advance_iop_forcing(const Spack& scm_dt, const Spack& ps_in, const uview_1d<const Spack>& u_in, const uview_1d<const Spack>& v_in, const uview_1d<const Spack>& t_in, const uview_1d<const Spack>& q_in, const uview_1d<const Spack>& t_phys_frc, const uview_1d<Spack>& u_update, const uview_1d<Spack>& v_update, const uview_1d<Spack>& t_update, const uview_1d<Spack>& q_update)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
