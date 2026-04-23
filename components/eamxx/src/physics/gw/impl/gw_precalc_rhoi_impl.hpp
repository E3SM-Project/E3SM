#ifndef GW_GWD_PRECALC_RHOI_IMPL_HPP
#define GW_GWD_PRECALC_RHOI_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_precalc_rhoi. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_precalc_rhoi(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const Real& dt,
  const Int& tend_level,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& t,
  const uview_2d<const Real>& gwut,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& nm,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& c,
  const uview_2d<const Real>& q,
  const uview_1d<const Real>& dse,
  // Outputs
  const uview_1d<Real>& egwdffi,
  const uview_2d<Real>& qtgw,
  const uview_1d<Real>& dttdf,
  const uview_1d<Real>& dttke,
  const uview_1d<Real>& ttgw)
{


  int nan_count = 0, inf_count = 0;
  int m_npgw = pgwv * 2 + 1;


  
  // rhoi_kludge: Recalculated rhoi to preserve answers.
  uview_1d<Real> rhoi_kludge, decomp_ca, decomp_cc, decomp_dnom, decomp_ze, q_nostride, qtgw_nostride;
  workspace.template take_many_contiguous_unsafe<7>(
    {"rhoi_kludge", "decomp_ca", "decomp_cc", "decomp_dnom", "decomp_ze", "q_nostride", "qtgw_nostride"},
    {&rhoi_kludge, &decomp_ca, &decomp_cc, &decomp_dnom, &decomp_ze, &q_nostride, &qtgw_nostride});

  rhoi_kludge(0) = pint(0) / (C::Rair.value * t(0));
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, 1, pver), [&] (const int k) {
      rhoi_kludge(k) = pint(k) * 2 / (C::Rair.value * (t(k) + t(k-1)));
    });
  rhoi_kludge(pver) = pint(pver) / (C::Rair.value * t(pver-1));


  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    for (int pg = 0; pg < m_npgw; pg++) {
      if (Kokkos::isnan(gwut(k,pg))) nan_count++;
      if (Kokkos::isinf(gwut(k,pg))) inf_count++;
    }
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gwd_precalc_rhoi] gwut before gw_ediff NaN=%d Inf=%d\n", nan_count, inf_count);
  }



  // Calculate effective diffusivity and LU decomposition for the
  // vertical diffusion solver.
  gw_ediff (team, workspace, pver, pgwv, init.kbotbg, init.ktop, tend_level, dt,
            gwut, ubm, nm, rhoi_kludge, pmid, rdpm, c,
            egwdffi, decomp_ca, decomp_cc, decomp_dnom, decomp_ze);


  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    for (int pg = 0; pg < m_npgw; pg++) {
      if (Kokkos::isnan(gwut(k,pg))) nan_count++;
      if (Kokkos::isinf(gwut(k,pg))) inf_count++;
    }
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gwd_precalc_rhoi] gwut after gw_ediff NaN=%d Inf=%d\n", nan_count, inf_count);
  }


  // Calculate tendency on each constituent.
  const int pcnst = q.extent(1);
  for (int m = 0; m < pcnst; ++m) {
    auto q_stride    = ekat::subview_1(q, m);
    auto qtgw_stride = ekat::subview_1(qtgw, m);
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
        q_nostride(k) = q_stride(k);
      });
    team.team_barrier();
    gw_diff_tend(team, workspace, pver, init.kbotbg, init.ktop, q_nostride, dt,
                 decomp_ca, decomp_cc, decomp_dnom, decomp_ze, qtgw_nostride);
    team.team_barrier();
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
        qtgw_stride(k) = qtgw_nostride(k);
      });
  }

  // Calculate tendency from diffusing dry static energy (dttdf).
  gw_diff_tend(team, workspace, pver, init.kbotbg, init.ktop, dse, dt, decomp_ca, decomp_cc, decomp_dnom, decomp_ze, dttdf);







  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    if (Kokkos::isnan(decomp_ze(k))) nan_count++;
    if (Kokkos::isinf(decomp_ze(k))) inf_count++;
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gwd_precalc_rhoi] decomp_ze NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    if (Kokkos::isnan(rhoi_kludge(k))) nan_count++;
    if (Kokkos::isinf(rhoi_kludge(k))) inf_count++;
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gwd_precalc_rhoi] rhoi_kludge NaN=%d Inf=%d\n", nan_count, inf_count);
  }

  nan_count = 0; inf_count = 0;
  for (int k = 0; k < pver; k++) {
    if (Kokkos::isnan(dttdf(k))) nan_count++;
    if (Kokkos::isinf(dttdf(k))) inf_count++;
  }
  if (nan_count > 0 || inf_count > 0) {
    printf("[gwd_precalc_rhoi] dttdf NaN=%d Inf=%d\n", nan_count, inf_count);
  }






  // Evaluate second temperature tendency term: Conversion of kinetic
  // energy into thermal.
  const int num_pgwv = 2*pgwv + 1;
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, (init.ktop+1), (init.kbotbg+1)), [&] (const int k) {
      for (int l = 0; l < num_pgwv; ++l) {
        dttke(k) += c(l) * gwut(k,l);
      }
    });

  team.team_barrier();
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      ttgw(k) = dttke(k) + dttdf(k);
    });

  workspace.template release_many_contiguous<7>(
    {&rhoi_kludge, &decomp_ca, &decomp_cc, &decomp_dnom, &decomp_ze, &q_nostride, &qtgw_nostride});
}

} // namespace gw
} // namespace scream

#endif
