#ifndef GW_GW_CM_SRC_IMPL_HPP
#define GW_GW_CM_SRC_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_cm_src. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_cm_src(
  // Inputs
  const MemberType& team,
  const GwCommonInit& init,
  const GwFrontInit& finit,
  const Int& pver,
  const Int& pgwv,
  const Int& kbot,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  const uview_1d<const Real>& frontgf,
  // Outputs
  Int& src_level,
  Int& tend_level,
  const uview_2d<Real>& tau,
  const uview_1d<Real>& ubm,
  const uview_1d<Real>& ubi,
  Real& xv,
  Real& yv,
  const uview_1d<Real>& c)
{
  //------------------------------------------------------------------------
  // Determine the source layer wind and unit vectors, then project winds.
  //------------------------------------------------------------------------
  gw_front_project_winds(team, pver, kbot, u, v, xv, yv, ubm, ubi);

  //-----------------------------------------------------------------------
  // Gravity wave sources
  //-----------------------------------------------------------------------
  gw_front_gw_sources(team, finit, pgwv, pver, kbot, frontgf, tau);

  src_level = kbot;
  tend_level = kbot;

  // Set phase speeds as reference speeds plus the wind speed at the source
  // level.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.cref.size()), [&] (const int l) {
      c(l) = init.cref(l) + std::abs(ubi(kbot+1));
    });
}

} // namespace gw
} // namespace scream

#endif
