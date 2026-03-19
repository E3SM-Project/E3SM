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

template<typename S, typename D>
void Functions<S,D>::gw_convect_init(
  // Inputs
  ekat::ParameterList& params,
  const Kokkos::View<Real***, Kokkos::HostSpace>& mfcc_in)
{
  s_convect_init.gw_convect_eff             = params.get<Real>("gw_convect_eff", s_convect_init.gw_convect_eff);
  s_convect_init.gw_convect_hcf             = params.get<Real>("gw_convect_hcf", s_convect_init.gw_convect_hcf);
  s_convect_init.gw_convect_hdepth_scale    = params.get<Real>("gw_convect_hdepth_scale", s_convect_init.gw_convect_hdepth_scale);
  s_convect_init.gw_convect_hdepth_min      = params.get<Real>("gw_convect_hdepth_min", s_convect_init.gw_convect_hdepth_min);
  s_convect_init.gw_convect_storm_speed_min = params.get<Real>("gw_convect_storm_speed_min", s_convect_init.gw_convect_storm_speed_min);
  s_convect_init.gw_convect_plev_src_wind   = params.get<Real>("gw_convect_plev_src_wind", s_convect_init.gw_convect_plev_src_wind);
  s_convect_init.use_gw_convect_old         = params.get<bool>("use_gw_convect_old", s_convect_init.use_gw_convect_old);
  s_convect_init.mfcc = view_3d<Real>("mfcc", mfcc_in.extent(0), mfcc_in.extent(1), mfcc_in.extent(2));
  Kokkos::deep_copy(s_convect_init.mfcc, mfcc_in);
}

} // namespace gw
} // namespace scream

#endif
