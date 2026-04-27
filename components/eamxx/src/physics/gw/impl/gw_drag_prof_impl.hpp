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

  int nan_count = 0, inf_count = 0;
  int m_npgw = pgwv * 2 + 1;


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

  //------------------------------------------------------------------------
  // Compute the stress profiles and diffusivities
  //------------------------------------------------------------------------
  gwd_compute_stress_profiles_and_diffusivities(
    team, workspace, init, pver, pgwv,
    src_level, ubi, c, rhoi,
    ni, kvtt, t, ti, piln, tau);

  // Tau projected in the four cardinal directions, for the momentum
  // conservation routine and for diagnostic output.
  if ( pgwv > 0) {
    gwd_project_tau(team, workspace, init, pver, pgwv, tend_level, tau, ubi, c, xv, yv, taucd);
  }

  //------------------------------------------------------------------------
  // Compute the tendencies from the stress divergence.
  //------------------------------------------------------------------------

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    if (Kokkos::isnan(ubm(k))) nan_count++;
    if (Kokkos::isinf(ubm(k))) inf_count++;
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gw_drag_prof] before compute ubm NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    for (int pg = 0; pg < m_npgw; pg++) {
      if (Kokkos::isnan(gwut(k,pg))) nan_count++;
      if (Kokkos::isinf(gwut(k,pg))) inf_count++;
    }
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gw_drag_prof] before compute gwut NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  gwd_compute_tendencies_from_stress_divergence(
    team, workspace, init, pver, pgwv, do_taper, dt, effgw, tend_level, max_level,
    lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw);

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    if (Kokkos::isnan(ubm(k))) nan_count++;
    if (Kokkos::isinf(ubm(k))) inf_count++;
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gw_drag_prof] after compute ubm NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    for (int pg = 0; pg < m_npgw; pg++) {
      if (Kokkos::isnan(gwut(k,pg))) nan_count++;
      if (Kokkos::isinf(gwut(k,pg))) inf_count++;
    }
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gw_drag_prof] after compute gwut NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  if (pgwv > 0) {
     // Precalculate rhoi for the following routines. We have rhoi, but
     // recalculate it here to preserve answers. (Ideal gas law.)
    gwd_precalc_rhoi(
      team, workspace, init, pver, pgwv, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse,
      egwdffi, qtgw, dttdf, dttke, ttgw);




    

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(egwdffi(k))) nan_count++;
      if (Kokkos::isinf(egwdffi(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[gw_drag_prof] egwdffi NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(qtgw(k,0))) nan_count++;
      if (Kokkos::isinf(qtgw(k,0))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[gw_drag_prof] qtgw NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(dttdf(k))) nan_count++;
      if (Kokkos::isinf(dttdf(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[gw_drag_prof] dttdf NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(dttke(k))) nan_count++;
      if (Kokkos::isinf(dttke(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[gw_drag_prof] dttke NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(ttgw(k))) nan_count++;
      if (Kokkos::isinf(ttgw(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[gw_drag_prof] ttgw NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    
    
    
    

  }
  else {
    const int pcnst = q.extent(1);
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, pver*pcnst), [&] (const int k_pcnst) {
        const int k = k_pcnst / pcnst;
        const int p = k_pcnst % pcnst;
        qtgw(k, p) = 0;
      });
  }
}

} // namespace gw
} // namespace scream

#endif
