#ifndef GW_GW_DRAG_PROF_IMPL_HPP
#define GW_GW_DRAG_PROF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_drag_prof. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_drag_prof(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const Int& src_level,
  const Int& max_level,
  const Int& tend_level,
  const bool& do_taper,
  const Real& dt,
  const Real& lat,
  const uview_1d<const Real>& t,
  const uview_1d<const Real>& ti,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& dpm,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& piln,
  const uview_1d<const Real>& rhoi,
  const uview_1d<const Real>& nm,
  const uview_1d<const Real>& ni,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& ubi,
  const Real& xv,
  const Real& yv,
  const Real& effgw,
  const uview_1d<const Real>& c,
  const uview_1d<const Real>& kvtt,
  const uview_2d<const Real>& q,
  const uview_1d<const Real>& dse,
  // Inputs/Outputs
  const uview_2d<Real>& tau,
  // Outputs
  const uview_1d<Real>& utgw,
  const uview_1d<Real>& vtgw,
  const uview_1d<Real>& ttgw,
  const uview_2d<Real>& qtgw,
  const uview_2d<Real>& taucd,
  const uview_1d<Real>& egwdffi,
  const uview_2d<Real>& gwut,
  const uview_1d<Real>& dttdf,
  const uview_1d<Real>& dttke)
{
  const bool dbg = (team.league_rank() == 0) && (team.team_rank() == 0);
  if (dbg) Kokkos::printf("[gw_drag_prof] enter pver=%d pgwv=%d src_level=%d max_level=%d tend_level=%d\n",
                          (int)pver, (int)pgwv, (int)src_level, (int)max_level, (int)tend_level);
  //------------------------------------------------------------------------
  // Initialize gravity wave drag tendencies to zero.
  //------------------------------------------------------------------------
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      utgw(k) = 0;
      vtgw(k) = 0;
      dttke(k) = 0;
      ttgw(k) = 0;
    });

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, (pver+1)*4), [&] (const int k_4) {
      const int k   = k_4 / 4;
      const int dir = k_4 % 4;
      taucd(k, dir) = 0;
    });

  const int num_pgwv = 2*pgwv + 1;
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver*num_pgwv), [&] (const int k_pgwv) {
      const int k = k_pgwv / num_pgwv;
      const int l = k_pgwv % num_pgwv;
      gwut(k, l) = 0;
    });

  team.team_barrier();
  if (dbg) Kokkos::printf("[gw_drag_prof] after zero-init parallel_fors\n");

  //------------------------------------------------------------------------
  // Compute the stress profiles and diffusivities
  //------------------------------------------------------------------------
  // DIAGNOSTIC: skip stress_prof entirely; just do a workspace take + release
  // to see if the hang is in workspace handling alone, or in stress_prof's
  // k-loop work.
  if (dbg) Kokkos::printf("[gw_drag_prof] DIAG: skipping stress_prof; testing workspace take/release only\n");
  {
    uview_1d<Real> ws_a, ws_b, ws_c, ws_d;
    workspace.template take_many_contiguous_unsafe<4>(
      {"diag_a", "diag_b", "diag_c", "diag_d"},
      {&ws_a, &ws_b, &ws_c, &ws_d});
    if (dbg) Kokkos::printf("[gw_drag_prof] DIAG: workspace take ok\n");
    workspace.template release_many_contiguous<4>({&ws_a, &ws_b, &ws_c, &ws_d});
    if (dbg) Kokkos::printf("[gw_drag_prof] DIAG: workspace release ok\n");
  }
  team.team_barrier();
  if (dbg) Kokkos::printf("[gw_drag_prof] after workspace diag test\n");

  // Tau projected in the four cardinal directions, for the momentum
  // conservation routine and for diagnostic output.
  if ( pgwv > 0) {
    gwd_project_tau(team, workspace, init, pver, pgwv, tend_level, tau, ubi, c, xv, yv, taucd);
    team.team_barrier();
    if (dbg) Kokkos::printf("[gw_drag_prof] after gwd_project_tau\n");
  } else {
    if (dbg) Kokkos::printf("[gw_drag_prof] skip gwd_project_tau (pgwv=0)\n");
  }

  //------------------------------------------------------------------------
  // Compute the tendencies from the stress divergence.
  //------------------------------------------------------------------------

  gwd_compute_tendencies_from_stress_divergence(
    team, workspace, init, pver, pgwv, do_taper, dt, effgw, tend_level, max_level,
    lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw);
  team.team_barrier();
  if (dbg) Kokkos::printf("[gw_drag_prof] after gwd_compute_tendencies_from_stress_divergence\n");

  if (pgwv > 0) {
     // Precalculate rhoi for the following routines. We have rhoi, but
     // recalculate it here to preserve answers. (Ideal gas law.)
    gwd_precalc_rhoi(
      team, workspace, init, pver, pgwv, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse,
      egwdffi, qtgw, dttdf, dttke, ttgw);
    team.team_barrier();
    if (dbg) Kokkos::printf("[gw_drag_prof] after gwd_precalc_rhoi\n");
  }
  else {
    const int pcnst = q.extent(1);
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, pver*pcnst), [&] (const int k_pcnst) {
        const int k = k_pcnst / pcnst;
        const int p = k_pcnst % pcnst;
        qtgw(k, p) = 0;
      });
    team.team_barrier();
    if (dbg) Kokkos::printf("[gw_drag_prof] after qtgw zero-out (pgwv=0 branch)\n");
  }
  if (dbg) Kokkos::printf("[gw_drag_prof] exit\n");
}

} // namespace gw
} // namespace scream

#endif
