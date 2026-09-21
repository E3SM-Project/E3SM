#ifndef GW_GW_BERES_SRC_IMPL_HPP
#define GW_GW_BERES_SRC_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_beres_src. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_beres_src(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const GwConvectInit& cinit,
  const Int& pver,
  const Int& pgwv,
  const Real& lat,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  const uview_1d<const Real>& netdt,
  const uview_1d<const Real>& zm,
  const Real& maxq0_conversion_factor,
  const Real& hdepth_scaling_factor,
  const Real& hdepth_min,
  const Real& storm_speed_min,
  const bool& use_gw_convect_old,
  // Outputs
  Int& src_level,
  Int& tend_level,
  const uview_2d<Real>& tau,
  const uview_1d<Real>& ubm,
  const uview_1d<Real>& ubi,
  Real& xv,
  Real& yv,
  const uview_1d<Real>& c,
  Real& hdepth,
  Real& maxq0_out)
{
  // Note all locals must be shared per team so that we have team-consistent
  // views of these values.
  auto local_storage = workspace.take("local_storage");

  // These alias shared workspace, so every write to them (here and in the
  // helper routines below) must come from a single thread.

  // Maximum heating rate.
  Real& maxq0 = local_storage(0);

  // Bottom/top heating range index.
  Int& mini = reinterpret_cast<Int&>(local_storage(1));
  Int& maxi = reinterpret_cast<Int&>(local_storage(2));

  // Mean wind in heating region.
  Real& uh = local_storage(3);

  // Min/max projected wind value in each column.
  Real& Umin = local_storage(4);
  Real& Umax = local_storage(5);

  // Speed of convective cells relative to storm.
  Int& storm_speed = reinterpret_cast<Int&>(local_storage(6));

  Kokkos::single(Kokkos::PerTeam(team), [&] {
    maxq0 = 0;
    mini = 0;
    maxi = 0;
    uh = 0;
    Umin = 0;
    Umax = 0;
    storm_speed = 0;
  });
  team.team_barrier();

  // note: the heating_altitude_max is probably not needed because there is
  // rarely any convective heating above this level and the performance impact
  // of skipping the iteration over higher levels is likely negilible.

  //----------------------------------------------------------------------
  // Initialize tau array
  //----------------------------------------------------------------------

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, tau.size()), [&] (const int k) {
      tau.data()[k] = 0;
    });

  team.team_barrier();

  // hdepth is a shared output (a view element in GWDrag::run_impl), so it
  // must have a single writer; gw_heating_depth likewise sets it from a
  // Kokkos::single.
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    hdepth = 0;
  });


  //------------------------------------------------------------------------
  // Determine source layer wind and unit vectors, then project winds.
  //------------------------------------------------------------------------

  gw_convect_project_winds(team, cinit, pver, u, v, xv, yv, ubm, ubi);

  //-----------------------------------------------------------------------
  // Calculate heating depth.
  //
  // Heating depth is defined as the first height range from the bottom in
  // which heating rate is continuously positive.
  //-----------------------------------------------------------------------
  gw_heating_depth(team, cinit, pver, maxq0_conversion_factor, hdepth_scaling_factor,
                   use_gw_convect_old, zm, netdt, mini, maxi, hdepth, maxq0_out, maxq0);

  //-----------------------------------------------------------------------
  // Taking ubm at assumed source level to be the storm speed,
  // find the cell speed where the storm speed is > storm_speed_min
  //-----------------------------------------------------------------------
  gw_storm_speed(team, init, cinit, pver, storm_speed_min, ubm, mini, maxi,
                 storm_speed, uh, Umin, Umax);

  //-----------------------------------------------------------------------
  // Gravity wave sources
  //-----------------------------------------------------------------------
  gw_convect_gw_sources(team, init, cinit, pgwv, pver, lat, hdepth_min, hdepth, mini, maxi, netdt, uh, storm_speed, maxq0, Umin, Umax, tau);

  // Output the source level.
  src_level = maxi;
  tend_level = maxi;

  // Set phase speeds; just use reference speeds.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.cref.size()), [&] (const int l) {
      c(l) = init.cref(l);
    });

  // mini/maxi/uh/maxq0/storm_speed alias this slot and were read just above
  team.team_barrier();

  workspace.release(local_storage);
}

} // namespace gw
} // namespace scream

#endif
