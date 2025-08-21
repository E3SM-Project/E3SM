#include "gw_test_data.hpp"

#include <ekat_math_utils.hpp>
#include <ekat_kokkos_types.hpp>
#include <ekat_assert.hpp>
#include <ekat_subview_utils.hpp>

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to gw fortran calls and vice versa
//

namespace scream {
namespace gw {

using GWF = Functions<Real, DefaultDevice>;
using GWC = typename GWF::C;

using ExeSpace   = typename GWF::KT::ExeSpace;
using MemberType = typename GWF::KT::MemberType;

using view1di_d = GWF::view_1d<Int>;
using view1dr_d = GWF::view_1d<Real>;
using view2dr_d = GWF::view_2d<Real>;
using view3dr_d = GWF::view_3d<Real>;

using WSM = typename GWF::WorkspaceManager;

extern "C" {

void gwd_compute_tendencies_from_stress_divergence_c(Int ncol, bool do_taper, Real dt, Real effgw, Int* tend_level, Real* lat, Real* dpm, Real* rdpm, Real* c, Real* ubm, Real* t, Real* nm, Real* xv, Real* yv, Real* tau, Real* gwut, Real* utgw, Real* vtgw);

void gw_init_c(Int pver_in, Int pgwv_in, Real dc_in, Real* cref_in, bool orographic_only, bool do_molec_diff_in, bool tau_0_ubc_in, Int nbot_molec_in, Int ktop_in, Int kbotbg_in, Real fcrit2_in, Real kwv_in, Real gravit_in, Real rair_in, Real* alpha_in);

void gw_prof_c(Int ncol, Real cpair, Real* t, Real* pmid, Real* pint, Real* rhoi, Real* ti, Real* nm, Real* ni);

void momentum_energy_conservation_c(Int ncol, Int* tend_level, Real dt, Real* taucd, Real* pint, Real* pdel, Real* u, Real* v, Real* dudt, Real* dvdt, Real* dsdt, Real* utgw, Real* vtgw, Real* ttgw);

void gwd_compute_stress_profiles_and_diffusivities_c(Int ncol, Int* src_level, Real* ubi, Real* c, Real* rhoi, Real* ni, Real* kvtt, Real* t, Real* ti, Real* piln, Real* tau);

void gwd_project_tau_c(Int ncol, Int* tend_level, Real* tau, Real* ubi, Real* c, Real* xv, Real* yv, Real* taucd);

void gwd_precalc_rhoi_c(Int pcnst, Int ncol, Real dt, Int* tend_level, Real* pmid, Real* pint, Real* t, Real* gwut, Real* ubm, Real* nm, Real* rdpm, Real* c, Real* q, Real* dse, Real* egwdffi, Real* qtgw, Real* dttdf, Real* dttke, Real* ttgw);

void gw_drag_prof_c(Int pcnst, Int ncol, Int* src_level, Int* tend_level, bool do_taper, Real dt, Real* lat, Real* t, Real* ti, Real* pmid, Real* pint, Real* dpm, Real* rdpm, Real* piln, Real* rhoi, Real* nm, Real* ni, Real* ubm, Real* ubi, Real* xv, Real* yv, Real effgw, Real* c, Real* kvtt, Real* q, Real* dse, Real* tau, Real* utgw, Real* vtgw, Real* ttgw, Real* qtgw, Real* taucd, Real* egwdffi, Real* gwut, Real* dttdf, Real* dttke);

void gw_front_init_c(Real taubgnd, Real frontgfc_in, Int kfront_in);

void gw_front_project_winds_c(Int ncol, Int kbot, Real* u, Real* v, Real* xv, Real* yv, Real* ubm, Real* ubi);

void gw_front_gw_sources_c(Int ncol, Int kbot, Real* frontgf, Real* tau);

void gw_cm_src_c(Int ncol, Int kbot, Real* u, Real* v, Real* frontgf, Int* src_level, Int* tend_level, Real* tau, Real* ubm, Real* ubi, Real* xv, Real* yv, Real* c);

void gw_convect_init_c(Int maxh, Int maxuh, Real plev_src_wind, Real* mfcc_in);

void gw_convect_project_winds_c(Int ncol, Real* u, Real* v, Real* xv, Real* yv, Real* ubm, Real* ubi);

void gw_heating_depth_c(Int ncol, Real maxq0_conversion_factor, Real hdepth_scaling_factor, bool use_gw_convect_old, Real* zm, Real* netdt, Int* mini, Int* maxi, Real* hdepth, Real* maxq0_out, Real* maxq0);

void gw_storm_speed_c(Int ncol, Real storm_speed_min, Real* ubm, Int* mini, Int* maxi, Int* storm_speed, Real* uh, Real* umin, Real* umax);

void gw_convect_gw_sources_c(Int ncol, Real* lat, Real hdepth_min, Real* hdepth, Int* mini, Int* maxi, Real* netdt, Real* uh, Int* storm_speed, Real* maxq0, Real* umin, Real* umax, Real* tau);

void gw_beres_src_c(Int ncol, Real* lat, Real* u, Real* v, Real* netdt, Real* zm, Int* src_level, Int* tend_level, Real* tau, Real* ubm, Real* ubi, Real* xv, Real* yv, Real* c, Real* hdepth, Real* maxq0_out, Real maxq0_conversion_factor, Real hdepth_scaling_factor, Real hdepth_min, Real storm_speed_min, bool use_gw_convect_old);

void gw_ediff_c(Int ncol, Int kbot, Int ktop, Int* tend_level, Real* gwut, Real* ubm, Real* nm, Real* rho, Real dt, Real gravit, Real* pmid, Real* rdpm, Real* c, Real* egwdffi, Real *decomp_ca, Real *decomp_cc, Real *decomp_dnom, Real *decomp_ze);

void gw_diff_tend_c(Int ncol, Int kbot, Int ktop, Real* q, Real dt, Real *decomp_ca, Real *decomp_cc, Real *decomp_dnom, Real *decomp_ze, Real* dq);

void gw_oro_init_c();

void gw_oro_src_c(Int ncol, Real* u, Real* v, Real* t, Real* sgh, Real* pmid, Real* pint, Real* dpm, Real* zm, Real* nm, Int* src_level, Int* tend_level, Real* tau, Real* ubm, Real* ubi, Real* xv, Real* yv, Real* c);

} // extern "C" : end _c decls

// Wrapper around gw_init
void gw_init(GwInit& init)
{
  gw_init_c(init.pver, init.pgwv, init.dc, init.cref, init.orographic_only, init.do_molec_diff, init.tau_0_ubc, init.nbot_molec, init.ktop, init.kbotbg, init.fcrit2, init.kwv, GWC::gravit, GWC::Rair, init.alpha);
}

// Wrapper around gw_init for cxx
void gw_init_cxx(GwInit& init)
{
  using uview_1d = typename GWF::uview_1d<Real>;
  GWF::gw_common_init(
    init.pver,
    init.pgwv,
    init.dc,
    uview_1d(init.cref, init.pgwv*2 + 1),
    init.orographic_only,
    init.do_molec_diff,
    init.tau_0_ubc,
    init.nbot_molec,
    init.ktop,
    init.kbotbg,
    init.fcrit2,
    init.kwv,
    uview_1d(init.alpha, init.pver + 1));
}

void gw_finalize_cxx(GwInit& init)
{
  GWF::gw_common_finalize();
}

void gwd_compute_tendencies_from_stress_divergence_f(GwdComputeTendenciesFromStressDivergenceData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>(); // This will shift array data + 1
  gwd_compute_tendencies_from_stress_divergence_c(d.ncol, d.do_taper, d.dt, d.effgw, d.tend_level, d.lat, d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.xv, d.yv, d.tau, d.gwut, d.utgw, d.vtgw);
  d.transition<ekat::TransposeDirection::f2c>(); // This will shift array data - 1
}

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d)
{
  gw_init_cxx(d.init);

  // create device views and copy
  std::vector<view1di_d> one_d_ints_in(1);
  std::vector<view1dr_d> one_d_reals_in(3);
  std::vector<view2dr_d> two_d_reals_in(8);
  std::vector<view3dr_d> three_d_reals_in(2);

  ekat::host_to_device({d.tend_level}, d.ncol, one_d_ints_in);
  ekat::host_to_device({d.lat, d.xv, d.yv}, d.ncol, one_d_reals_in);
  ekat::host_to_device({d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.utgw, d.vtgw},
                       std::vector<int>(8, d.ncol),
                       std::vector<int>{    // dim2 sizes
                         d.init.pver,       // dpm
                         d.init.pver,       // rdpm
                         2*d.init.pgwv + 1, // c
                         d.init.pver,       // ubm
                         d.init.pver,       // t
                         d.init.pver,       // nm
                         d.init.pver,       // utgw
                         d.init.pver},      // vtgw
                       two_d_reals_in);
  ekat::host_to_device({d.tau, d.gwut},
                       std::vector<int>(2, d.ncol),
                       std::vector<int>{2*d.init.pgwv + 1, d.init.pver},
                       std::vector<int>{d.init.pver + 1, 2*d.init.pgwv + 1},
                       three_d_reals_in);

  const auto tend_level = one_d_ints_in[0];

  const auto lat        = one_d_reals_in[0];
  const auto xv         = one_d_reals_in[1];
  const auto yv         = one_d_reals_in[2];

  const auto dpm        = two_d_reals_in[0];
  const auto rdpm       = two_d_reals_in[1];
  const auto c          = two_d_reals_in[2];
  const auto ubm        = two_d_reals_in[3];
  const auto t          = two_d_reals_in[4];
  const auto nm         = two_d_reals_in[5];
  const auto utgw       = two_d_reals_in[6];
  const auto vtgw       = two_d_reals_in[7];

  const auto tau        = three_d_reals_in[0];
  const auto gwut       = three_d_reals_in[1];

  // Find max tend_level
  int max_level = 0;
  Kokkos::parallel_reduce("find max level", d.ncol, KOKKOS_LAMBDA(const int i, int& lmax) {
    if (tend_level(i) > lmax) {
      lmax = tend_level(i);
    }
  }, Kokkos::Max<int>(max_level));

  auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.ncol, d.init.pver);

  WSM wsm((d.init.pver + 1) * (2*d.init.pgwv + 1), 1, policy);
  GWF::GwCommonInit init_cp = GWF::s_common_init;

  // unpack init because we do not want the lambda to capture it
  const int pver = d.init.pver;
  const int pgwv = d.init.pgwv;
  const bool do_taper = d.do_taper;
  const Real dt = d.dt;
  const Real effgw = d.effgw;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int col = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto dpm_c  = ekat::subview(dpm, col);
    const auto rdpm_c = ekat::subview(rdpm, col);
    const auto c_c    = ekat::subview(c, col);
    const auto ubm_c  = ekat::subview(ubm, col);
    const auto t_c    = ekat::subview(t, col);
    const auto nm_c   = ekat::subview(nm, col);
    const auto tau_c  = ekat::subview(tau, col);
    const auto utgw_c = ekat::subview(utgw, col);
    const auto vtgw_c = ekat::subview(vtgw, col);
    const auto gwut_c = ekat::subview(gwut, col);

    GWF::gwd_compute_tendencies_from_stress_divergence(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver, pgwv, do_taper, dt, effgw,
      tend_level(col),
      max_level,
      lat(col),
      dpm_c,
      rdpm_c,
      c_c,
      ubm_c,
      t_c,
      nm_c,
      xv(col),
      yv(col),
      tau_c,
      gwut_c,
      utgw_c,
      vtgw_c
    );
  });

  // Get outputs back
  std::vector<view2dr_d> two_d_reals_out = {utgw, vtgw};
  std::vector<view3dr_d> three_d_reals_out = {tau, gwut};
  ekat::device_to_host({d.utgw, d.vtgw}, d.ncol, d.init.pver, two_d_reals_out);
  ekat::device_to_host({d.tau, d.gwut},
                       std::vector<int>(2, d.ncol),
                       std::vector<int>{2*d.init.pgwv + 1, d.init.pver},
                       std::vector<int>{d.init.pver + 1, 2*d.init.pgwv + 1},
                       three_d_reals_out);

  gw_finalize_cxx(d.init);
}

void gw_prof(GwProfData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_prof_c(d.ncol, d.cpair, d.t, d.pmid, d.pint, d.rhoi, d.ti, d.nm, d.ni);
  d.transition<ekat::TransposeDirection::f2c>();
}

void momentum_energy_conservation_f(MomentumEnergyConservationData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  momentum_energy_conservation_c(d.ncol, d.tend_level, d.dt, d.taucd, d.pint, d.pdel, d.u, d.v, d.dudt, d.dvdt, d.dsdt, d.utgw, d.vtgw, d.ttgw);
  d.transition<ekat::TransposeDirection::f2c>();
}

void momentum_energy_conservation(MomentumEnergyConservationData& d)
{
  gw_init_cxx(d.init);

  // create device views and copy
  std::vector<view1di_d> one_d_ints_in(1);
  std::vector<view2dr_d> two_d_reals_in(10);
  std::vector<view3dr_d> three_d_reals_in(1);

  ekat::host_to_device({d.tend_level}, d.ncol, one_d_ints_in);
  std::vector<int> dim2_sizes(10, d.init.pver);
  dim2_sizes[0] += 1; // pint
  ekat::host_to_device({d.pint, d.pdel, d.u, d.v, d.dudt, d.dvdt, d.dsdt, d.utgw, d.vtgw, d.ttgw},
                       std::vector<int>(10, d.ncol),
                       dim2_sizes,
                       two_d_reals_in);
  ekat::host_to_device({d.taucd},
                       d.ncol, d.init.pver + 1, 4,
                       three_d_reals_in);

  const auto tend_level = one_d_ints_in[0];

  const auto pint = two_d_reals_in[0];
  const auto pdel = two_d_reals_in[1];
  const auto u    = two_d_reals_in[2];
  const auto v    = two_d_reals_in[3];
  const auto dudt = two_d_reals_in[4];
  const auto dvdt = two_d_reals_in[5];
  const auto dsdt = two_d_reals_in[6];
  const auto utgw = two_d_reals_in[7];
  const auto vtgw = two_d_reals_in[8];
  const auto ttgw = two_d_reals_in[9];

  const auto taucd = three_d_reals_in[0];

  auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.ncol, d.init.pver);

  WSM wsm((d.init.pver + 1) * (2*d.init.pgwv + 1), 1, policy);
  GWF::GwCommonInit init_cp = GWF::s_common_init;

  // unpack init because we do not want the lambda to capture it
  const int pver = d.init.pver;
  const Real dt = d.dt;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int col = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto pint_c  = ekat::subview(pint, col);
    const auto pdel_c  = ekat::subview(pdel, col);
    const auto u_c     = ekat::subview(u, col);
    const auto v_c     = ekat::subview(v, col);
    const auto dudt_c  = ekat::subview(dudt, col);
    const auto dvdt_c  = ekat::subview(dvdt, col);
    const auto dsdt_c  = ekat::subview(dsdt, col);
    const auto utgw_c  = ekat::subview(utgw, col);
    const auto vtgw_c  = ekat::subview(vtgw, col);
    const auto ttgw_c  = ekat::subview(ttgw, col);
    const auto taucd_c = ekat::subview(taucd, col);

    GWF::momentum_energy_conservation(
      team,
      pver,
      tend_level(col),
      dt,
      taucd_c,
      pint_c,
      pdel_c,
      u_c,
      v_c,
      dudt_c,
      dvdt_c,
      dsdt_c,
      utgw_c,
      vtgw_c,
      ttgw_c
    );
  });

  // Get outputs back
  std::vector<view2dr_d> two_d_reals_out = {dudt, dvdt, dsdt, utgw, vtgw, ttgw};
  ekat::device_to_host({d.dudt, d.dvdt, d.dsdt, d.utgw, d.vtgw, d.ttgw},
                       d.ncol, d.init.pver, two_d_reals_out);

  gw_finalize_cxx(d.init);
}

void gwd_compute_stress_profiles_and_diffusivities_f(GwdComputeStressProfilesAndDiffusivitiesData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gwd_compute_stress_profiles_and_diffusivities_c(d.ncol, d.src_level, d.ubi, d.c, d.rhoi, d.ni, d.kvtt, d.t, d.ti, d.piln, d.tau);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gwd_compute_stress_profiles_and_diffusivities(GwdComputeStressProfilesAndDiffusivitiesData& d)
{
  gw_init_cxx(d.init);

  // create device views and copy
  std::vector<view1di_d> one_d_ints_in(1);
  std::vector<view2dr_d> two_d_reals_in(8);
  std::vector<view3dr_d> three_d_reals_in(1);

  ekat::host_to_device({d.src_level}, d.ncol, one_d_ints_in);
  ekat::host_to_device({d.ubi, d.c, d.rhoi, d.ni, d.kvtt, d.t, d.ti, d.piln},
                       std::vector<int>(8, d.ncol),
                       std::vector<int>{    // dim2 sizes
                         d.init.pver + 1,   // ubi
                         2*d.init.pgwv + 1, // c
                         d.init.pver + 1,   // rhoi
                         d.init.pver + 1,   // ni
                         d.init.pver + 1,   // kvtt
                         d.init.pver,       // t
                         d.init.pver + 1,   // ti
                         d.init.pver + 1},  // piln
                       two_d_reals_in);
  ekat::host_to_device({d.tau}, d.ncol, 2*d.init.pgwv + 1, d.init.pver + 1, three_d_reals_in);

  const auto src_level = one_d_ints_in[0];

  const auto ubi  = two_d_reals_in[0];
  const auto c    = two_d_reals_in[1];
  const auto rhoi = two_d_reals_in[2];
  const auto ni   = two_d_reals_in[3];
  const auto kvtt = two_d_reals_in[4];
  const auto t    = two_d_reals_in[5];
  const auto ti   = two_d_reals_in[6];
  const auto piln = two_d_reals_in[7];

  const auto tau        = three_d_reals_in[0];

  auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.ncol, d.init.pver);

  WSM wsm((d.init.pver + 1) * (2*d.init.pgwv + 1), 4, policy);
  GWF::GwCommonInit init_cp = GWF::s_common_init;

  // unpack init because we do not want the lambda to capture it
  const int pver = d.init.pver;
  const int pgwv = d.init.pgwv;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int col = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto ubi_c  = ekat::subview(ubi, col);
    const auto c_c    = ekat::subview(c, col);
    const auto rhoi_c = ekat::subview(rhoi, col);
    const auto ni_c   = ekat::subview(ni, col);
    const auto kvtt_c = ekat::subview(kvtt, col);
    const auto t_c    = ekat::subview(t, col);
    const auto ti_c   = ekat::subview(ti, col);
    const auto piln_c = ekat::subview(piln, col);
    const auto tau_c  = ekat::subview(tau, col);

    GWF::gwd_compute_stress_profiles_and_diffusivities(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver, pgwv,
      src_level(col),
      ubi_c,
      c_c,
      rhoi_c,
      ni_c,
      kvtt_c,
      t_c,
      ti_c,
      piln_c,
      tau_c);
  });

  // Get outputs back
  std::vector<view3dr_d> three_d_reals_out = {tau};
  ekat::device_to_host({d.tau}, d.ncol, 2*d.init.pgwv + 1, d.init.pver + 1, three_d_reals_out);

  gw_finalize_cxx(d.init);
}

void gwd_project_tau(GwdProjectTauData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gwd_project_tau_c(d.ncol, d.tend_level, d.tau, d.ubi, d.c, d.xv, d.yv, d.taucd);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gwd_precalc_rhoi(GwdPrecalcRhoiData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gwd_precalc_rhoi_c(d.pcnst, d.ncol, d.dt, d.tend_level, d.pmid, d.pint, d.t, d.gwut, d.ubm, d.nm, d.rdpm, d.c, d.q, d.dse, d.egwdffi, d.qtgw, d.dttdf, d.dttke, d.ttgw);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_drag_prof(GwDragProfData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_drag_prof_c(d.pcnst, d.ncol, d.src_level, d.tend_level, d.do_taper, d.dt, d.lat, d.t, d.ti, d.pmid, d.pint, d.dpm, d.rdpm, d.piln, d.rhoi, d.nm, d.ni, d.ubm, d.ubi, d.xv, d.yv, d.effgw, d.c, d.kvtt, d.q, d.dse, d.tau, d.utgw, d.vtgw, d.ttgw, d.qtgw, d.taucd, d.egwdffi, d.gwut, d.dttdf, d.dttke);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_front_init(GwFrontInitData& d)
{
  gw_init(d.init);
  gw_front_init_c(d.taubgnd, d.frontgfc_in, d.kfront_in);
}

void gw_front_project_winds(GwFrontProjectWindsData& d)
{
  gw_front_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_front_project_winds_c(d.ncol, d.kbot, d.u, d.v, d.xv, d.yv, d.ubm, d.ubi);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_front_gw_sources(GwFrontGwSourcesData& d)
{
  gw_front_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_front_gw_sources_c(d.ncol, d.kbot, d.frontgf, d.tau);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_cm_src(GwCmSrcData& d)
{
  gw_front_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_cm_src_c(d.ncol, d.kbot, d.u, d.v, d.frontgf, d.src_level, d.tend_level, d.tau, d.ubm, d.ubi, d.xv, d.yv, d.c);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_convect_init(GwConvectInitData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_convect_init_c(d.maxh, d.maxuh, d.plev_src_wind, d.mfcc_in);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_convect_project_winds(GwConvectProjectWindsData& d)
{
  gw_convect_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_convect_project_winds_c(d.ncol, d.u, d.v, d.xv, d.yv, d.ubm, d.ubi);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_heating_depth(GwHeatingDepthData& d)
{
  gw_convect_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_heating_depth_c(d.ncol, d.maxq0_conversion_factor, d.hdepth_scaling_factor, d.use_gw_convect_old, d.zm, d.netdt, d.mini, d.maxi, d.hdepth, d.maxq0_out, d.maxq0);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_storm_speed(GwStormSpeedData& d)
{
  gw_convect_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_storm_speed_c(d.ncol, d.storm_speed_min, d.ubm, d.mini, d.maxi, d.storm_speed, d.uh, d.umin, d.umax);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_convect_gw_sources(GwConvectGwSourcesData& d)
{
  gw_convect_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_convect_gw_sources_c(d.ncol, d.lat, d.hdepth_min, d.hdepth, d.mini, d.maxi, d.netdt, d.uh, d.storm_speed, d.maxq0, d.umin, d.umax, d.tau);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_beres_src(GwBeresSrcData& d)
{
  gw_convect_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_beres_src_c(d.ncol, d.lat, d.u, d.v, d.netdt, d.zm, d.src_level, d.tend_level, d.tau, d.ubm, d.ubi, d.xv, d.yv, d.c, d.hdepth, d.maxq0_out, d.maxq0_conversion_factor, d.hdepth_scaling_factor, d.hdepth_min, d.storm_speed_min, d.use_gw_convect_old);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_ediff(GwEdiffData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_ediff_c(d.ncol, d.kbot, d.ktop, d.tend_level, d.gwut, d.ubm, d.nm, d.rho, d.dt, GWC::gravit, d.pmid, d.rdpm, d.c, d.egwdffi, d.decomp_ca, d.decomp_cc, d.decomp_dnom, d.decomp_ze);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_diff_tend(GwDiffTendData& d)
{
  gw_init(d.init);
  d.transition<ekat::TransposeDirection::c2f>();
  gw_diff_tend_c(d.ncol, d.kbot, d.ktop, d.q, d.dt, d.decomp_ca, d.decomp_cc, d.decomp_dnom, d.decomp_ze, d.dq);
  d.transition<ekat::TransposeDirection::f2c>();
}

void gw_oro_src(GwOroSrcData& d)
{
  gw_init(d.init);
  gw_oro_init_c();
  d.transition<ekat::TransposeDirection::c2f>();
  gw_oro_src_c(d.ncol, d.u, d.v, d.t, d.sgh, d.pmid, d.pint, d.dpm, d.zm, d.nm, d.src_level, d.tend_level, d.tau, d.ubm, d.ubi, d.xv, d.yv, d.c);
  d.transition<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace gw
} // namespace scream
