#ifndef GW_GW_CONVECT_GW_SOURCES_IMPL_HPP
#define GW_GW_CONVECT_GW_SOURCES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_convect_gw_sources. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_convect_gw_sources(
  // Inputs
  const MemberType& team,
  const GwCommonInit& init,
  const GwConvectInit& cinit,
  const Int& pgwv,
  const Int& pver,
  const Real& lat,
  const Real& hdepth_min,
  const Real& hdepth,
  const Int& mini,
  const Int& maxi,
  const uview_1d<const Real>& netdt,
  const Real& uh,
  const Int& storm_speed,
  const Real& maxq0,
  const Real& umin,
  const Real& umax,
  // Outputs
  const uview_2d<Real>& tau)
{
  const Int num_pgwv = 2*pgwv + 1;

  //---------------------------------------------------------------------
  // Look up spectrum only if depth >= 2.5 km, else set tau0 = 0.
  //---------------------------------------------------------------------
  if ((hdepth >= hdepth_min) && (std::abs(lat) < (C::Pi/2))) {
    //------------------------------------------------------------------
    // Look up the spectrum using depth and uh.
    //------------------------------------------------------------------

    // Shift spectrum so that it is relative to the ground.
    const Int shift = -static_cast<Int>(std::round(storm_speed/init.dc)) % num_pgwv;

    // Adjust for critical level filtering.
    const Int Umini = ekat::impl::max(static_cast<Int>(std::round(umin/init.dc)), -pgwv) + pgwv;
    const Int Umaxi = ekat::impl::min(static_cast<Int>(std::round(umax/init.dc)), pgwv) + pgwv;

    const Int hdepth_i = static_cast<Int>(std::round(hdepth)) - 1;
    const Int uh_i = static_cast<Int>(std::round(uh)) + cinit.maxuh;
    EKAT_KERNEL_ASSERT_MSG(hdepth_i >= 0 && hdepth_i < cinit.maxh, "Bad hdepth");
    EKAT_KERNEL_ASSERT_MSG(uh_i >= 0 && uh_i < 2*cinit.maxuh + 1, "Bad uh");
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, num_pgwv), [&] (const int l) {
        if (Umaxi > Umini && (l >= Umini && l <= Umaxi)) {
          tau(l, maxi+1) = 0;
        }
        else {
          Int shifted_l = (l + shift + num_pgwv) % num_pgwv;
          tau(l, maxi+1) = cinit.mfcc(hdepth_i, uh_i, shifted_l) * maxq0 * maxq0 / GWC::tau_avg_length;
        }
      });
  }  // depth >= 2.5
}

} // namespace gw
} // namespace scream

#endif
