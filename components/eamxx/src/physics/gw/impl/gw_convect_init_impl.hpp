#ifndef GW_GW_CONVECT_INIT_IMPL_HPP
#define GW_GW_CONVECT_INIT_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_convect_init. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::gw_convect_init(
  // Inputs
  const Real& plev_src_wind,
  const uview_3d<const Real>& mfcc_in)
{
  // Just set k_src_wind to pver. We don't have access to pref_edge
  s_convect_init.k_src_wind = s_common_init.pver - 1;

  // First dimension is maxh.
  s_convect_init.maxh   = mfcc_in.extent(0);

  // Second dimension is -maxuh:maxuh (size 2*maxuh+1).
  s_convect_init.maxuh = (mfcc_in.extent(1) - 1) / 2;

  s_convect_init.mfcc = view_3d<Real>("mfcc", mfcc_in.extent(0), mfcc_in.extent(1), mfcc_in.extent(2));
  Kokkos::deep_copy(s_convect_init.mfcc, mfcc_in);
}

} // namespace gw
} // namespace scream

#endif
