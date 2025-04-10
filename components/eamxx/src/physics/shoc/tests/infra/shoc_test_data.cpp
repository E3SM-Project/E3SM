#include "shoc_test_data.hpp"

#include "shoc_data.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include "share/util/eamxx_deep_copy.hpp"

#include <random>

using scream::Real;
using scream::Int;

namespace scream {
namespace shoc {

//
// Glue functions to call from host with the Data struct
//
// We are provisionally following P3 here in case SHOC uses global data.
//

void shoc_grid(ShocGridData& d)
{
  shoc_grid_host(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt, d.dz_zi, d.rho_zt);
}

void shoc_diag_obklen(ShocDiagObklenData& d)
{
  shoc_diag_obklen_host(d.shcol, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc, d.thl_sfc, d.cldliq_sfc, d.qv_sfc, d.ustar, d.kbfs, d.obklen);
}

void update_host_dse(UpdateHostDseData& d)
{
  update_host_dse_host(d.shcol, d.nlev, d.thlm, d.shoc_ql, d.inv_exner, d.zt_grid, d.phis, d.host_dse);
}

void shoc_energy_fixer(ShocEnergyFixerData& d)
{
  shoc_energy_fixer_host(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.zt_grid, d.zi_grid, d.se_b, d.ke_b, d.wv_b, d.wl_b, d.se_a, d.ke_a, d.wv_a, d.wl_a, d.wthl_sfc, d.wqw_sfc, d.rho_zt, d.tke, d.pint, d.host_dse);
}

void shoc_energy_integrals(ShocEnergyIntegralsData& d)
{
  shoc_energy_integrals_host(d.shcol, d.nlev, d.host_dse, d.pdel, d.rtm, d.rcm, d.u_wind, d.v_wind, d.se_int, d.ke_int, d.wv_int, d.wl_int);
}

void calc_shoc_vertflux(CalcShocVertfluxData& d)
{
  calc_shoc_vertflux_host(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar, d.vertflux);
}

void calc_shoc_varorcovar(CalcShocVarorcovarData& d)
{
  calc_shoc_varorcovar_host(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi, d.dz_zi, d.invar1, d.invar2, d.varorcovar);
}

void compute_tmpi(ComputeTmpiData& d)
{
  compute_tmpi_host(d.nlevi, d.shcol, d.dtime, d.rho_zi, d.dz_zi, d.tmpi);
}

void dp_inverse(DpInverseData& d)
{
  dp_inverse_host(d.nlev, d.shcol, d.rho_zt, d.dz_zt, d.rdp_zt);
}

void integ_column_stability(IntegColumnStabilityData& d)
{
  integ_column_stability_host(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
}

void check_tke(CheckTkeData& d)
{
  check_tke_host(d.shcol, d.nlev, d.tke);
}

void shoc_tke(ShocTkeData& d)
{
  shoc_tke_host(d.shcol, d.nlev, d.nlevi, d.dtime, d.wthv_sec, d.shoc_mix, d.dz_zi, d.dz_zt, d.pres, d.tabs, d.u_wind, d.v_wind, d.brunt, d.zt_grid, d.zi_grid, d.pblh, d.tke, d.tk, d.tkh, d.isotropy);
}

void compute_shr_prod(ComputeShrProdData& d)
{
  compute_shr_prod_host(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind, d.v_wind, d.sterm);
}

void isotropic_ts(IsotropicTsData& d)
{
  isotropic_ts_host(d.nlev, d.shcol, d.brunt_int, d.tke, d.a_diss, d.brunt, d.isotropy);
}

void adv_sgs_tke(AdvSgsTkeData& d)
{
  adv_sgs_tke_host(d.nlev, d.shcol, d.dtime, d.shoc_mix, d.wthv_sec, d.sterm_zt, d.tk, d.tke, d.a_diss);
}

void eddy_diffusivities(EddyDiffusivitiesData& d)
{
  eddy_diffusivities_host(d.nlev, d.shcol, d.pblh, d.zt_grid, d.tabs, d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
}

void shoc_length(ShocLengthData& d)
{
  shoc_length_host(d.shcol, d.nlev, d.nlevi, d.host_dx, d.host_dy, d.zt_grid, d.zi_grid, d.dz_zt, d.tke, d.thv, d.brunt, d.shoc_mix);
}

void compute_brunt_shoc_length(ComputeBruntShocLengthData& d)
{
  compute_brunt_shoc_length_host(d.nlev, d.nlevi, d.shcol, d.dz_zt, d.thv, d.thv_zi, d.brunt);
}

void compute_l_inf_shoc_length(ComputeLInfShocLengthData& d)
{
  compute_l_inf_shoc_length_host(d.nlev, d.shcol, d.zt_grid, d.dz_zt, d.tke, d.l_inf);
}

void compute_shoc_mix_shoc_length(ComputeShocMixShocLengthData& d)
{
  compute_shoc_mix_shoc_length_host(d.nlev, d.shcol, d.tke, d.brunt, d.zt_grid, d.l_inf, d.shoc_mix);
}

void check_length_scale_shoc_length(CheckLengthScaleShocLengthData& d)
{
  check_length_scale_shoc_length_host(d.nlev, d.shcol, d.host_dx, d.host_dy, d.shoc_mix);
}

void clipping_diag_third_shoc_moments(ClippingDiagThirdShocMomentsData& d)
{
  clipping_diag_third_shoc_moments_host(d.nlevi, d.shcol, d.w_sec_zi, d.w3);
}

void diag_second_moments_srf(DiagSecondMomentsSrfData& d)
{
  shoc_diag_second_moments_srf_host(d.shcol, d.wthl_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar);
}

void linear_interp(LinearInterpData& d)
{
  linear_interp_host(d.x1, d.x2, d.y1, d.y2, d.km1, d.km2, d.ncol, d.minthresh);
}

void diag_third_shoc_moments(DiagThirdShocMomentsData& d)
{
  diag_third_shoc_moments_host(d.shcol, d.nlev, d.nlevi, d.w_sec, d.thl_sec, d.wthl_sec, d.isotropy, d.brunt, d.thetal, d.tke, d.dz_zt, d.dz_zi, d.zt_grid, d.zi_grid, d.w3);
}

void compute_diag_third_shoc_moment(ComputeDiagThirdShocMomentData& d)
{
  compute_diag_third_shoc_moment_host(d.shcol, d.nlev, d.nlevi, d.w_sec, d.thl_sec, d.wthl_sec, d.tke, d.dz_zt, d.dz_zi, d.isotropy_zi, d.brunt_zi, d.w_sec_zi, d.thetal_zi, d.w3);
}

void shoc_assumed_pdf(ShocAssumedPdfData& d)
{
  shoc_assumed_pdf_host(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.w_field, d.thl_sec, d.qw_sec, d.wthl_sec, d.w_sec, d.wqw_sec, d.qwthl_sec, d.w3, d.pres, d.zt_grid, d.zi_grid, d.shoc_cldfrac, d.shoc_ql, d.wqls, d.wthv_sec, d.shoc_ql2);
}

void shoc_assumed_pdf_tilde_to_real(ShocAssumedPdfTildeToRealData& d)
{
  shoc_assumed_pdf_tilde_to_real_host(d.w_first, d.sqrtw2, &d.w1);
}

void shoc_assumed_pdf_vv_parameters(ShocAssumedPdfVvParametersData& d)
{
  shoc_assumed_pdf_vv_parameters_host(d.w_first, d.w_sec, d.w3var, d.w_tol_sqd, &d.skew_w, &d.w1_1, &d.w1_2, &d.w2_1, &d.w2_2, &d.a);
}

void shoc_assumed_pdf_thl_parameters(ShocAssumedPdfThlParametersData& d)
{
  shoc_assumed_pdf_thl_parameters_host(d.wthlsec, d.sqrtw2, d.sqrtthl, d.thlsec, d.thl_first, d.w1_1, d.w1_2, d.skew_w, d.a, d.thl_tol, d.w_thresh, &d.thl1_1, &d.thl1_2, &d.thl2_1, &d.thl2_2, &d.sqrtthl2_1, &d.sqrtthl2_2);
}

void shoc_assumed_pdf_qw_parameters(ShocAssumedPdfQwParametersData& d)
{
  shoc_assumed_pdf_qw_parameters_host(d.wqwsec, d.sqrtw2, d.skew_w, d.sqrtqt, d.qwsec, d.w1_2, d.w1_1, d.qw_first, d.a, d.rt_tol, d.w_thresh, &d.qw1_1, &d.qw1_2, &d.qw2_1, &d.qw2_2, &d.sqrtqw2_1, &d.sqrtqw2_2);
}

void shoc_assumed_pdf_inplume_correlations(ShocAssumedPdfInplumeCorrelationsData& d)
{
  shoc_assumed_pdf_inplume_correlations_host(d.sqrtqw2_1, d.sqrtthl2_1, d.a, d.sqrtqw2_2, d.sqrtthl2_2, d.qwthlsec, d.qw1_1, d.qw_first, d.thl1_1, d.thl_first, d.qw1_2, d.thl1_2, &d.r_qwthl_1);
}

void shoc_assumed_pdf_compute_temperature(ShocAssumedPdfComputeTemperatureData& d)
{
  shoc_assumed_pdf_compute_temperature_host(d.thl1, d.pval, &d.tl1);
}

void shoc_assumed_pdf_compute_qs(ShocAssumedPdfComputeQsData& d)
{
  shoc_assumed_pdf_compute_qs_host(d.tl1_1, d.tl1_2, d.pval, &d.qs1, &d.beta1, &d.qs2, &d.beta2);
}

void shoc_assumed_pdf_compute_s(ShocAssumedPdfComputeSData& d)
{
  shoc_assumed_pdf_compute_s_host(d.qw1, d.qs1, d.beta, d.pval, d.thl2, d.qw2, d.sqrtthl2, d.sqrtqw2, d.r_qwthl, &d.s, &d.std_s, &d.qn, &d.c);
}

void shoc_assumed_pdf_compute_sgs_liquid(ShocAssumedPdfComputeSgsLiquidData& d)
{
  shoc_assumed_pdf_compute_sgs_liquid_host(d.a, d.ql1, d.ql2, &d.shoc_ql);
}

void shoc_assumed_pdf_compute_cloud_liquid_variance(ShocAssumedPdfComputeCloudLiquidVarianceData& d)
{
  shoc_assumed_pdf_compute_cloud_liquid_variance_host(d.a, d.s1, d.ql1, d.c1, d.std_s1, d.s2, d.ql2, d.c2, d.std_s2, d.shoc_ql, &d.shoc_ql2);
}

void shoc_assumed_pdf_compute_liquid_water_flux(ShocAssumedPdfComputeLiquidWaterFluxData& d)
{
  shoc_assumed_pdf_compute_liquid_water_flux_host(d.a, d.w1_1, d.w_first, d.ql1, d.w1_2, d.ql2, &d.wqls);
}

void shoc_assumed_pdf_compute_buoyancy_flux(ShocAssumedPdfComputeBuoyancyFluxData& d)
{
  shoc_assumed_pdf_compute_buoyancy_flux_host(d.wthlsec, d.wqwsec, d.pval, d.wqls, &d.wthv_sec);
}

void diag_second_moments_ubycond(DiagSecondMomentsUbycondData& d)
{
  shoc_diag_second_moments_ubycond_host(d.shcol, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec);
}

void pblintd_init_pot(PblintdInitPotData& d)
{
  shoc_pblintd_init_pot_host(d.shcol, d.nlev, d.thl, d.ql, d.q, d.thv);
}

void pblintd_cldcheck(PblintdCldcheckData& d)
{
  shoc_pblintd_cldcheck_host(d.shcol, d.nlev, d.nlevi, d.zi, d.cldn, d.pblh);
}

void diag_second_moments_lbycond(DiagSecondMomentsLbycondData& d)
{
  diag_second_moments_lbycond_host(d.shcol, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar, d.wthl_sec, d.wqw_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.thl_sec, d.qw_sec, d.qwthl_sec);
}

void diag_second_moments(DiagSecondMomentsData& d)
{
  diag_second_moments_host(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh, d.tk,
                           d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec,
                           d.vw_sec, d.wtke_sec, d.w_sec);
}

void diag_second_shoc_moments(DiagSecondShocMomentsData& d)
{
  diag_second_shoc_moments_host(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh,
                                d.tk, d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.thl_sec, d.qw_sec,
                                d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.w_sec);
}

void compute_shoc_vapor(ComputeShocVaporData& d)
{
  compute_shoc_vapor_host(d.shcol, d.nlev, d.qw, d.ql, d.qv);
}

void update_prognostics_implicit(UpdatePrognosticsImplicitData& d)
{
  update_prognostics_implicit_host(d.shcol, d.nlev, d.nlevi, d.num_tracer, d.dtime,
                                   d.dz_zt, d.dz_zi, d.rho_zt, d.zt_grid, d.zi_grid,
                                   d.tk, d.tkh, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc,
                                   d.wtracer_sfc, d.thetal, d.qw, d.tracer, d.tke, d.u_wind, d.v_wind);
}

void shoc_main(ShocMainData& d)
{
  const int npbl = shoc_init_host(d.nlev, d.pref_mid, d.nbot_shoc, d.ntop_shoc);
  d.elapsed_s = shoc_main_host(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, npbl, d.host_dx, d.host_dy, d.thv, d.zt_grid, d.zi_grid,
              d.pres, d.presi, d.pdel, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.wtracer_sfc,
              d.num_qtracers, d.w_field, d.inv_exner, d.phis, d.host_dse, d.tke, d.thetal, d.qw,
              d.u_wind, d.v_wind, d.qtracers, d.wthv_sec, d.tkh, d.tk, d.shoc_ql, d.shoc_cldfrac, d.pblh,
              d.shoc_mix, d.isotropy, d.w_sec, d.thl_sec, d.qw_sec, d.qwthl_sec, d.wthl_sec, d.wqw_sec,
              d.wtke_sec, d.uw_sec, d.vw_sec, d.w3, d.wqls_sec, d.brunt, d.shoc_ql2);
}

void pblintd_height(PblintdHeightData& d)
{
  pblintd_height_host(d.shcol, d.nlev, d.npbl, d.z, d.u, d.v, d.ustar, d.thv, d.thv_ref, d.pblh, d.rino, d.check);
}

void vd_shoc_decomp_and_solve(VdShocDecompandSolveData& d)
{
  vd_shoc_decomp_and_solve_host(d.shcol, d.nlev, d.nlevi, d.n_rhs, d.dtime, d.kv_term, d.tmpi, d.rdp_zt, d.flux, d.var);
}

void pblintd_surf_temp(PblintdSurfTempData& d)
{
  pblintd_surf_temp_host(d.shcol, d.nlev, d.nlevi, d.z, d.ustar, d.obklen, d.kbfs, d.thv, d.tlv, d.pblh, d.check, d.rino);
}

void pblintd_check_pblh(PblintdCheckPblhData& d)
{
  pblintd_check_pblh_host(d.shcol, d.nlev, d.nlevi, d.nlev/*npbl*/, d.z, d.ustar, d.check, d.pblh);
}

void pblintd(PblintdData& d)
{
  pblintd_host(d.shcol, d.nlev, d.nlevi, d.npbl, d.z, d.zi, d.thl, d.ql, d.q, d.u, d.v, d.ustar, d.obklen, d.kbfs, d.cldn, d.pblh);
}

void compute_shoc_temperature(ComputeShocTempData& d)
{
  compute_shoc_temperature_host(d.shcol, d.nlev, d.thetal, d.ql, d.inv_exner, d.tabs);
}

// end _c impls

//
// _host function definitions. These expect data in C layout
//

void calc_shoc_varorcovar_host(Int shcol, Int nlev, Int nlevi, Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
                            Real *invar1, Real *invar2, Real *varorcovar)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 6;

  std::vector<view_2d> temp_d(num_arrays);
  std::vector<int> dim1_sizes(num_arrays, shcol);
  std::vector<int> dim2_sizes        = {nlevi,       nlevi,  nlevi, nlev,   nlev,   nlevi};
  std::vector<const Real*> ptr_array = {isotropy_zi, tkh_zi, dz_zi, invar1, invar2, varorcovar};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    isotropy_zi_d  (temp_d[0]),
    tkh_zi_d    (temp_d[1]),
    dz_zi_d     (temp_d[2]),
    invar1_d    (temp_d[3]),
    invar2_d    (temp_d[4]),
    varorcovar_d(temp_d[5]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto oisotropy_zi_d = ekat::subview(isotropy_zi_d, i);
    const auto otkh_zi_d      = ekat::subview(tkh_zi_d, i);
    const auto odz_zi_d       = ekat::subview(dz_zi_d, i);
    const auto oinvar1_d      = ekat::subview(invar1_d, i);
    const auto oinvar2_d      = ekat::subview(invar2_d, i);
    const auto ovarorcovar_d  = ekat::subview(varorcovar_d, i);

    SHF::calc_shoc_varorcovar(team, nlev, tunefac, oisotropy_zi_d, otkh_zi_d, odz_zi_d,
                              oinvar1_d, oinvar2_d, ovarorcovar_d);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {varorcovar_d};
  ekat::device_to_host({varorcovar}, shcol, nlevi, inout_views);
}

void calc_shoc_vertflux_host(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
        Real *dz_zi, Real *invar, Real *vertflux)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 4;

  std::vector<view_2d> temp_d(num_arrays);
  std::vector<Int> dim1_sizes(num_arrays, shcol);
  std::vector<Int> dim2_sizes        = {nlevi,  nlevi, nlev,  nlevi};
  std::vector<const Real*> ptr_array = {tkh_zi, dz_zi, invar, vertflux};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    tkh_zi_d  (temp_d[0]),
    dz_zi_d   (temp_d[1]),
    invar_d   (temp_d[2]),
    vertflux_d(temp_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto otkh_zi_d   = ekat::subview(tkh_zi_d, i);
    const auto odz_zi_d    = ekat::subview(dz_zi_d, i);
    const auto oinvar_d    = ekat::subview(invar_d, i);
    const auto overtflux_d = ekat::subview(vertflux_d, i);

    SHF::calc_shoc_vertflux(team, nlev, otkh_zi_d, odz_zi_d, oinvar_d, overtflux_d);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {vertflux_d};
  ekat::device_to_host({vertflux}, shcol, nlevi, inout_views);
}

void shoc_diag_second_moments_srf_host(Int shcol, Real* wthl_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using view_1d    = typename SHOC::view_1d<Scalar>;

  std::vector<view_1d> temp_d(3);
  ScreamDeepCopy::copy_to_device({wthl_sfc, uw_sfc, vw_sfc}, shcol, temp_d);

  // inputs
  view_1d
    wthl_sfc_d (temp_d[0]),
    uw_sfc_d   (temp_d[1]),
    vw_sfc_d   (temp_d[2]);

  // outputs
  view_1d ustar2_d("ustar2", shcol),
          wstar_d ("wstar", shcol);

  Kokkos::parallel_for("parallel_moments_srf", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar wthl_sfc_s{wthl_sfc_d(i)};
     Scalar uw_sfc_s{uw_sfc_d(i)};
     Scalar vw_sfc_s{vw_sfc_d(i)};

     Scalar ustar2_s{0};
     Scalar wstar_s{0};

     SHOC::shoc_diag_second_moments_srf(wthl_sfc_s, uw_sfc_s, vw_sfc_s, ustar2_s, wstar_s);

     ustar2_d(i) = ustar2_s;
     wstar_d(i)  = wstar_s;
   });

  std::vector<view_1d> inout_views = {ustar2_d, wstar_d};
  ScreamDeepCopy::copy_to_host({ustar2, wstar}, shcol, inout_views);
}

void shoc_diag_second_moments_ubycond_host(Int shcol, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec,
      Real* wtke_sec)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using view_1d    = typename SHOC::view_1d<Scalar>;

  view_1d thl_sec_d  ("thl_sec"  ,shcol),
          qw_sec_d   ("qw_sec"   ,shcol),
          qwthl_sec_d("qwthl_sec",shcol),
          wthl_sec_d ("wthl_sec" ,shcol),
          wqw_sec_d  ("wqw_sec"  ,shcol),
          uw_sec_d   ("uw_sec"   ,shcol),
          vw_sec_d   ("vw_sec"   ,shcol),
          wtke_sec_d ("wtke_sec" ,shcol);

  Kokkos::parallel_for("parallel_moments_ubycond", shcol, KOKKOS_LAMBDA (const int& i) {

    Scalar thl_sec_s{0.};
    Scalar qw_sec_s{0.};
    Scalar wthl_sec_s{0.};
    Scalar wqw_sec_s{0.};
    Scalar qwthl_sec_s{0.};
    Scalar uw_sec_s{0.};
    Scalar vw_sec_s{0.};
    Scalar wtke_sec_s{0.};

    SHOC::shoc_diag_second_moments_ubycond(thl_sec_s, qw_sec_s, wthl_sec_s, wqw_sec_s, qwthl_sec_s, uw_sec_s, vw_sec_s, wtke_sec_s);

    thl_sec_d(i)   = thl_sec_s;
    qw_sec_d(i)    = qw_sec_s;
    wthl_sec_d(i)  = wthl_sec_s;
    wqw_sec_d(i)   = wqw_sec_s;
    qwthl_sec_d(i) = qwthl_sec_s;
    uw_sec_d(i)    = uw_sec_s;
    vw_sec_d(i)    = vw_sec_s;
    wtke_sec_d(i)  = wtke_sec_s;

  });

  std::vector<view_1d> host_views = {thl_sec_d, qw_sec_d, qwthl_sec_d, wthl_sec_d, wqw_sec_d, uw_sec_d, vw_sec_d, wtke_sec_d};

  ScreamDeepCopy::copy_to_host({thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec}, shcol, host_views);
}

void update_host_dse_host(Int shcol, Int nlev, Real* thlm, Real* shoc_ql, Real* inv_exner, Real* zt_grid,
                       Real* phis, Real* host_dse)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_1d> temp_1d_d(1);
  std::vector<view_2d> temp_2d_d(5);
  std::vector<const Real*> ptr_array = {thlm,  shoc_ql, inv_exner, zt_grid, host_dse};

  // Sync to device
  ScreamDeepCopy::copy_to_device({phis}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d);

  view_1d phis_d(temp_1d_d[0]);

  view_2d
    thlm_d     (temp_2d_d[0]),
    shoc_ql_d  (temp_2d_d[1]),
    inv_exner_d(temp_2d_d[2]),
    zt_grid_d  (temp_2d_d[3]),
    host_dse_d (temp_2d_d[4]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar phis_s{phis_d(i)};
    const auto thlm_s      = ekat::subview(thlm_d, i);
    const auto shoc_ql_s   = ekat::subview(shoc_ql_d, i);
    const auto inv_exner_s = ekat::subview(inv_exner_d, i);
    const auto zt_grid_s   = ekat::subview(zt_grid_d, i);
    const auto host_dse_s  = ekat::subview(host_dse_d, i);

    SHF::update_host_dse(team, nlev, thlm_s, shoc_ql_s, inv_exner_s, zt_grid_s, phis_s, host_dse_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {host_dse_d};
  ekat::device_to_host({host_dse}, shcol, nlev, inout_views);
}

void compute_diag_third_shoc_moment_host(Int shcol, Int nlev, Int nlevi, Real* w_sec,
                                      Real* thl_sec, Real* wthl_sec, Real* tke,
                                      Real* dz_zt, Real* dz_zi, Real* isotropy_zi,
                                      Real* brunt_zi, Real* w_sec_zi, Real* thetal_zi,
                                      Real* w3)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_2d> temp_d(11);
  std::vector<Int> dim1_sizes(11, shcol);
  std::vector<Int> dim2_sizes     = {nlev,        nlevi,
                                     nlevi,       nlev,
                                     nlev,        nlevi,
                                     nlevi,       nlevi,
                                     nlevi,       nlevi,
                                     nlevi};
  std::vector<const Real*> ptr_array = {w_sec,       thl_sec,
                                        wthl_sec,    tke,
                                        dz_zt,       dz_zi,
                                        isotropy_zi, brunt_zi,
                                        w_sec_zi,    thetal_zi,
                                        w3};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    w_sec_d      (temp_d[0]),
    thl_sec_d    (temp_d[1]),
    wthl_sec_d   (temp_d[2]),
    tke_d        (temp_d[3]),
    dz_zt_d      (temp_d[4]),
    dz_zi_d      (temp_d[5]),
    isotropy_zi_d(temp_d[6]),
    brunt_zi_d   (temp_d[7]),
    w_sec_zi_d   (temp_d[8]),
    thetal_zi_d  (temp_d[9]),
    w3_d         (temp_d[10]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto w_sec_s       = ekat::subview(w_sec_d, i);
    const auto thl_sec_s     = ekat::subview(thl_sec_d, i);
    const auto wthl_sec_s    = ekat::subview(wthl_sec_d, i);
    const auto tke_s         = ekat::subview(tke_d, i);
    const auto dz_zt_s       = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s       = ekat::subview(dz_zi_d, i);
    const auto isotropy_zi_s = ekat::subview(isotropy_zi_d, i);
    const auto brunt_zi_s    = ekat::subview(brunt_zi_d, i);
    const auto w_sec_zi_s    = ekat::subview(w_sec_zi_d, i);
    const auto thetal_zi_s   = ekat::subview(thetal_zi_d, i);
    const auto w3_s          = ekat::subview(w3_d, i);

    // Hardcode runtime options for F90 testing
    const Real c_diag_3rd_mom = 7.0;
    SHF::compute_diag_third_shoc_moment(team, nlev, nlevi, c_diag_3rd_mom, w_sec_s, thl_sec_s,
                                        wthl_sec_s, tke_s, dz_zt_s, dz_zi_s, isotropy_zi_s,
                                        brunt_zi_s, w_sec_zi_s, thetal_zi_s, w3_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {w3_d};
  ekat::device_to_host({w3}, shcol, nlevi, inout_views);
}

void shoc_pblintd_init_pot_host(Int shcol, Int nlev, Real *thl, Real* ql, Real* q,
                             Real *thv)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  static constexpr Int num_arrays = 3;

  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device({thl, ql, q}, shcol, nlev, temp_d);

  view_2d thl_d(temp_d[0]),
          ql_d (temp_d[1]),
          q_d  (temp_d[2]);

  view_2d thv_d("thv", shcol, nlev);

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto thl_1d = ekat::subview(thl_d, i);
    const auto ql_1d  = ekat::subview(ql_d, i);
    const auto q_1d   = ekat::subview(q_d, i);
    const auto thv_1d = ekat::subview(thv_d, i);

    SHOC::shoc_pblintd_init_pot(team, nlev, thl_1d, ql_1d, q_1d, thv_1d);
  });

  std::vector<view_2d> inout_views = {thv_d};
  ekat::device_to_host({thv}, shcol, nlev, inout_views);
}

void compute_shoc_mix_shoc_length_host(Int nlev, Int shcol, Real* tke, Real* brunt,
                                    Real* zt_grid, Real* l_inf, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_1d> temp_1d_d(1);
  std::vector<view_2d> temp_2d_d(4);
  std::vector<const Real*> ptr_array = {tke,   brunt, zt_grid, shoc_mix};

  // Sync to device
  ScreamDeepCopy::copy_to_device({l_inf}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d);

  view_1d
    l_inf_d  (temp_1d_d[0]);

  view_2d
    tke_d     (temp_2d_d[0]),
    brunt_d   (temp_2d_d[1]),
    zt_grid_d (temp_2d_d[2]),
    shoc_mix_d  (temp_2d_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar l_inf_s {l_inf_d(i)};

    const auto tke_s      = ekat::subview(tke_d, i);
    const auto brunt_s    = ekat::subview(brunt_d, i);
    const auto zt_grid_s  = ekat::subview(zt_grid_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    const Real length_fac = 0.5;
    SHF::compute_shoc_mix_shoc_length(team, nlev, length_fac, tke_s, brunt_s, zt_grid_s, l_inf_s,
                                      shoc_mix_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {shoc_mix_d};
  ekat::device_to_host({shoc_mix}, shcol, nlev, inout_views);
}

void check_tke_host(Int shcol, Int nlev, Real* tke)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  std::vector<view_2d> temp_2d_d(1);

  // Sync to device
  ekat::host_to_device({tke}, shcol, nlev, temp_2d_d);

  view_2d
    tke_d(temp_2d_d[0]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto tke_s   = ekat::subview(tke_d, i);

    SHOC::check_tke(team, nlev, tke_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tke_d};
  ekat::device_to_host({tke}, shcol, nlev, inout_views);
}

void linear_interp_host(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_2d> temp_2d_d(3);
  std::vector<Int> dim1_sizes(3, ncol);
  std::vector<Int> dim2_sizes     = {km1,  km2,  km1};
  std::vector<const Real*> ptr_array = {x1,   x2,   y1};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);

  view_2d
    x1_d(temp_2d_d[0]),
    x2_d(temp_2d_d[1]),
    y1_d(temp_2d_d[2]),
    y2_d("y2_d", ncol, km2);

  const Int nk_pack = ekat::npack<Spack>(km1);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto x1_s  = ekat::subview(x1_d, i);
    const auto x2_s  = ekat::subview(x2_d, i);
    const auto y1_s  = ekat::subview(y1_d, i);
    const auto y2_s  = ekat::subview(y2_d, i);

    SHF::linear_interp(team, x1_s, x2_s, y1_s, y2_s, km1, km2, minthresh);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {y2_d};
  ekat::device_to_host({y2}, ncol, km2, inout_views);
}

void clipping_diag_third_shoc_moments_host(Int nlevi, Int shcol, Real *w_sec_zi,
                                        Real *w3)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  // Sync to device
  std::vector<view_2d> temp_d(2);
  ekat::host_to_device({w_sec_zi, w3}, shcol, nlevi, temp_d);

  view_2d
    w_sec_zi_d(temp_d[0]),
    w3_d      (temp_d[1]);

  const Int nk_pack = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto w_sec_zi_s = ekat::subview(w_sec_zi_d, i);
    const auto w3_s       = ekat::subview(w3_d, i);

    SHF::clipping_diag_third_shoc_moments(team, nlevi, w_sec_zi_s, w3_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {w3_d};
  ekat::device_to_host({w3}, shcol, nlevi, inout_views);
}

void shoc_energy_integrals_host(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_2d> temp_d(6);
  std::vector<const Real*> ptr_array = {host_dse, pdel, rtm,   rcm,   u_wind, v_wind};

  // Sync to device
  ekat::host_to_device(ptr_array, shcol, nlev, temp_d);

  // inputs
  view_2d
    host_dse_d(temp_d[0]),
    pdel_d    (temp_d[1]),
    rtm_d     (temp_d[2]),
    rcm_d     (temp_d[3]),
    u_wind_d  (temp_d[4]),
    v_wind_d  (temp_d[5]);

  // outputs
  view_1d
    se_int_d("se_int", shcol),
    ke_int_d("ke_int", shcol),
    wv_int_d("wv_int", shcol),
    wl_int_d("wl_int", shcol);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto host_dse_s = ekat::subview(host_dse_d, i);
    const auto pdel_s     = ekat::subview(pdel_d, i);
    const auto rtm_s      = ekat::subview(rtm_d, i);
    const auto rcm_s      = ekat::subview(rcm_d, i);
    const auto u_wind_s   = ekat::subview(u_wind_d, i);
    const auto v_wind_s   = ekat::subview(v_wind_d, i);

    Scalar se_int_s{0};
    Scalar ke_int_s{0};
    Scalar wv_int_s{0};
    Scalar wl_int_s{0};

    SHF::shoc_energy_integrals(team, nlev, host_dse_s, pdel_s, rtm_s, rcm_s, u_wind_s, v_wind_s,
                               se_int_s, ke_int_s, wv_int_s, wl_int_s);

    se_int_d(i) = se_int_s;
    ke_int_d(i) = ke_int_s;
    wv_int_d(i) = wv_int_s;
    wl_int_d(i) = wl_int_s;
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {se_int_d, ke_int_d, wv_int_d, wl_int_d};
  ScreamDeepCopy::copy_to_host({se_int,ke_int,wv_int,wl_int}, shcol, inout_views);
}

void diag_second_moments_lbycond_host(Int shcol, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar,
     Real* wthl_sec, Real* wqw_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* thl_sec, Real* qw_sec, Real* qwthl_sec)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using view_1d    = typename SHOC::view_1d<Scalar>;

  std::vector<view_1d> lbycond_d(6);
  ScreamDeepCopy::copy_to_device({wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, ustar2, wstar}, shcol, lbycond_d);

  // inputs
  view_1d wthl_d  (lbycond_d[0]),
          wqw_d   (lbycond_d[1]),
          uw_d    (lbycond_d[2]),
          vw_d    (lbycond_d[3]),
          ustar2_d(lbycond_d[4]),
          wstar_d (lbycond_d[5]);

  // outputs
  view_1d wthlo_d ("wthl",  shcol),
          wqwo_d  ("wqw",   shcol),
          uwo_d   ("uw",    shcol),
          vwo_d   ("vw",    shcol),
          wtkeo_d ("wtke",  shcol),
          thlo_d  ("thl",   shcol),
          qwo_d   ("qw",    shcol),
          qwthlo_d("qwthl", shcol);

  Kokkos::parallel_for("parallel_moments_lbycond", shcol, KOKKOS_LAMBDA (const int& i) {

    Scalar wthl_s{wthl_d(i)},
           wqw_s{wqw_d(i)},
           uw_s{uw_d(i)},
           vw_s{vw_d(i)},
           ustar2_s{ustar2_d(i)},
           wstar_s{wstar_d(i)};

    Scalar wthlo_s{0.},
           wqwo_s{0.},
           uwo_s{0.},
           vwo_s{0.},
           wtkeo_s{0.},
           thlo_s{0.},
           qwo_s{0.},
           qwthlo_s{0.};

    SHOC::shoc_diag_second_moments_lbycond(wthl_s, wqw_s, uw_s, vw_s, ustar2_s, wstar_s,
                                          wthlo_s, wqwo_s, uwo_s, vwo_s, wtkeo_s, thlo_s, qwo_s, qwthlo_s);

    wthlo_d  (i) = wthlo_s;
    wqwo_d   (i) = wqwo_s;
    uwo_d    (i) = uwo_s;
    vwo_d    (i) = vwo_s;
    wtkeo_d  (i) = wtkeo_s;
    thlo_d   (i) = thlo_s;
    qwo_d    (i) = qwo_s;
    qwthlo_d (i) = qwthlo_s;
  });

  std::vector<view_1d> host_views = {wthlo_d, wqwo_d, uwo_d, vwo_d, wtkeo_d, thlo_d, qwo_d, qwthlo_d};
  ScreamDeepCopy::copy_to_host({wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec, thl_sec, qw_sec, qwthl_sec}, shcol, host_views);
}

void diag_second_moments_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind,
          Real* tke, Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix,
          Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec,
          Real* wtke_sec, Real* w_sec)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  std::vector<Int> dim1_array(20, shcol);
  std::vector<Int> dim2_array = {nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev, nlev,
                                 nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi};

  std::vector<view_2d> temp_2d(20);
  std::vector<const Real*> ptr_array = {thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, zt_grid, shoc_mix,
                      thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, dz_zi, zi_grid};

  ekat::host_to_device(ptr_array, dim1_array, dim2_array, temp_2d);

  view_2d
    thetal_2d   (temp_2d[0]),
    qw_2d       (temp_2d[1]),
    u_wind_2d   (temp_2d[2]),
    v_wind_2d   (temp_2d[3]),
    tke_2d      (temp_2d[4]),
    isotropy_2d (temp_2d[5]),
    tkh_2d      (temp_2d[6]),
    tk_2d       (temp_2d[7]),
    zt_grid_2d  (temp_2d[8]),
    shoc_mix_2d (temp_2d[9]),
    thl_sec_2d  (temp_2d[10]),
    qw_sec_2d   (temp_2d[11]),
    wthl_sec_2d (temp_2d[12]),
    wqw_sec_2d  (temp_2d[13]),
    qwthl_sec_2d(temp_2d[14]),
    uw_sec_2d   (temp_2d[15]),
    vw_sec_2d   (temp_2d[16]),
    wtke_sec_2d (temp_2d[17]),
    dz_zi_2d    (temp_2d[18]),
    zi_grid_2d  (temp_2d[19]);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  view_2d w_sec_2d("w_sec", shcol, nlev_packs),
          isotropy_zi_2d("isotropy_zi", shcol, nlevi_packs),
          tkh_zi_2d("tkh_zi", shcol, nlevi_packs),
          tk_zi_2d("tk_zi", shcol, nlevi_packs);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    // Hardcode runtime options for F90
    const Real thl2tune = 1.0;
    const Real qw2tune = 1.0;
    const Real qwthl2tune = 1.0;
    const Real w2tune = 1.0;

    const auto thetal_1d      = ekat::subview(thetal_2d, i);
    const auto qw_1d          = ekat::subview(qw_2d, i);
    const auto u_wind_1d      = ekat::subview(u_wind_2d, i);
    const auto v_wind_1d      = ekat::subview(v_wind_2d, i);
    const auto tke_1d         = ekat::subview(tke_2d, i);
    const auto isotropy_1d    = ekat::subview(isotropy_2d, i);
    const auto tkh_1d         = ekat::subview(tkh_2d, i);
    const auto tk_1d          = ekat::subview(tk_2d, i);
    const auto dz_zi_1d       = ekat::subview(dz_zi_2d, i);
    const auto zt_grid_1d     = ekat::subview(zt_grid_2d, i);
    const auto zi_grid_1d     = ekat::subview(zi_grid_2d, i);
    const auto shoc_mix_1d    = ekat::subview(shoc_mix_2d, i);
    const auto thl_sec_1d     = ekat::subview(thl_sec_2d, i);
    const auto qw_sec_1d      = ekat::subview(qw_sec_2d, i);
    const auto wthl_sec_1d    = ekat::subview(wthl_sec_2d, i);
    const auto wqw_sec_1d     = ekat::subview(wqw_sec_2d, i);
    const auto qwthl_sec_1d   = ekat::subview(qwthl_sec_2d, i);
    const auto uw_sec_1d      = ekat::subview(uw_sec_2d, i);
    const auto vw_sec_1d      = ekat::subview(vw_sec_2d, i);
    const auto wtke_sec_1d    = ekat::subview(wtke_sec_2d, i);
    const auto w_sec_1d       = ekat::subview(w_sec_2d, i);
    const auto isotropy_zi_1d = ekat::subview(isotropy_zi_2d, i);
    const auto tkh_zi_1d      = ekat::subview(tkh_zi_2d, i);
    const auto tk_zi_1d       = ekat::subview(tk_zi_2d, i);

    SHOC::diag_second_moments(team, nlev, nlevi, thl2tune, qw2tune, qwthl2tune, w2tune,
                     thetal_1d, qw_1d, u_wind_1d, v_wind_1d, tke_1d, isotropy_1d, tkh_1d, tk_1d,
                     dz_zi_1d, zt_grid_1d, zi_grid_1d, shoc_mix_1d, isotropy_zi_1d, tkh_zi_1d, tk_zi_1d,
                     thl_sec_1d, qw_sec_1d, wthl_sec_1d, wqw_sec_1d,
                     qwthl_sec_1d, uw_sec_1d, vw_sec_1d, wtke_sec_1d, w_sec_1d);


  });

  std::vector<Int> dim1(9, shcol);
  std::vector<Int> dim2 = {nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlev };
  std::vector<view_2d> host_views = {thl_sec_2d, qw_sec_2d, wthl_sec_2d, wqw_sec_2d, qwthl_sec_2d, uw_sec_2d, vw_sec_2d, wtke_sec_2d, w_sec_2d};
  ekat::device_to_host({thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec}, dim1, dim2, host_views);
}

void diag_second_shoc_moments_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke,
                                Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix,
                                Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* thl_sec, Real* qw_sec, Real* wthl_sec,
                                Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec)
{
  using SHOC          = Functions<Real, DefaultDevice>;
  using Scalar        = typename SHOC::Scalar;
  using Spack         = typename SHOC::Spack;
  using view_2d       = typename SHOC::view_2d<Spack>;
  using KT            = typename SHOC::KT;
  using ExeSpace      = typename KT::ExeSpace;
  using MemberType    = typename SHOC::MemberType;
  using view_1d = typename SHOC::view_1d<Scalar>;

  std::vector<view_1d> temp_1d(4);
  ScreamDeepCopy::copy_to_device({wthl_sfc, wqw_sfc, uw_sfc, vw_sfc}, shcol, temp_1d);

  view_1d wthl_1d  (temp_1d[0]),
          wqw_1d   (temp_1d[1]),
          uw_1d    (temp_1d[2]),
          vw_1d    (temp_1d[3]);

  std::vector<Int> dim1_array(20, shcol);
  std::vector<Int> dim2_array = {nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev,  nlev, nlev,
                                 nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi};
  std::vector<view_2d> temp_2d(20);
  std::vector<const Real*> ptr_array = {thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, zt_grid, shoc_mix,
                      thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, dz_zi, zi_grid};

  ekat::host_to_device(ptr_array, dim1_array, dim2_array, temp_2d);

  view_2d
    thetal_2d   (temp_2d[0]),
    qw_2d       (temp_2d[1]),
    u_wind_2d   (temp_2d[2]),
    v_wind_2d   (temp_2d[3]),
    tke_2d      (temp_2d[4]),
    isotropy_2d (temp_2d[5]),
    tkh_2d      (temp_2d[6]),
    tk_2d       (temp_2d[7]),
    zt_grid_2d  (temp_2d[8]),
    shoc_mix_2d (temp_2d[9]),
    thl_sec_2d  (temp_2d[10]),
    qw_sec_2d   (temp_2d[11]),
    wthl_sec_2d (temp_2d[12]),
    wqw_sec_2d  (temp_2d[13]),
    qwthl_sec_2d(temp_2d[14]),
    uw_sec_2d   (temp_2d[15]),
    vw_sec_2d   (temp_2d[16]),
    wtke_sec_2d (temp_2d[17]),
    dz_zi_2d    (temp_2d[18]),
    zi_grid_2d  (temp_2d[19]);

  view_2d w_sec_2d("w_sec", shcol, nlev);

  view_1d ustar2_1d("ustar2", shcol),
                wstar_1d("wstar", shcol);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 3, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    // Hardcode runtime options for F90
    const Real thl2tune = 1.0;
    const Real qw2tune = 1.0;
    const Real qwthl2tune = 1.0;
    const Real w2tune = 1.0;

    auto workspace = workspace_mgr.get_workspace(team);

    const auto thetal_1d      = ekat::subview(thetal_2d, i);
    const auto qw_1d          = ekat::subview(qw_2d, i);
    const auto u_wind_1d      = ekat::subview(u_wind_2d, i);
    const auto v_wind_1d      = ekat::subview(v_wind_2d, i);
    const auto tke_1d         = ekat::subview(tke_2d, i);
    const auto isotropy_1d    = ekat::subview(isotropy_2d, i);
    const auto tkh_1d         = ekat::subview(tkh_2d, i);
    const auto tk_1d          = ekat::subview(tk_2d, i);
    const auto dz_zi_1d       = ekat::subview(dz_zi_2d, i);
    const auto zt_grid_1d     = ekat::subview(zt_grid_2d, i);
    const auto zi_grid_1d     = ekat::subview(zi_grid_2d, i);
    const auto shoc_mix_1d    = ekat::subview(shoc_mix_2d, i);
    const auto thl_sec_1d     = ekat::subview(thl_sec_2d, i);
    const auto qw_sec_1d      = ekat::subview(qw_sec_2d, i);
    const auto wthl_sec_1d    = ekat::subview(wthl_sec_2d, i);
    const auto wqw_sec_1d     = ekat::subview(wqw_sec_2d, i);
    const auto qwthl_sec_1d   = ekat::subview(qwthl_sec_2d, i);
    const auto uw_sec_1d      = ekat::subview(uw_sec_2d, i);
    const auto vw_sec_1d      = ekat::subview(vw_sec_2d, i);
    const auto wtke_sec_1d    = ekat::subview(wtke_sec_2d, i);
    const auto w_sec_1d       = ekat::subview(w_sec_2d, i);

    Scalar wthl_s   = wthl_1d(i);
    Scalar wqw_s    = wqw_1d(i);
    Scalar uw_s     = uw_1d(i);
    Scalar vw_s     = vw_1d(i);
    Scalar ustar2_s = ustar2_1d(i);
    Scalar wstar_s  = wstar_1d(i);

    SHOC::diag_second_shoc_moments(team, nlev, nlevi,
       thl2tune, qw2tune, qwthl2tune, w2tune,
       thetal_1d, qw_1d, u_wind_1d, v_wind_1d, tke_1d, isotropy_1d, tkh_1d, tk_1d, dz_zi_1d, zt_grid_1d, zi_grid_1d, shoc_mix_1d,
       wthl_s, wqw_s, uw_s, vw_s, ustar2_s, wstar_s,
       workspace, thl_sec_1d, qw_sec_1d, wthl_sec_1d, wqw_sec_1d, qwthl_sec_1d,
       uw_sec_1d, vw_sec_1d, wtke_sec_1d, w_sec_1d);
  });

  std::vector<Int> dim1(9, shcol);
  std::vector<Int> dim2 = {nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlev };
  std::vector<view_2d> host_2d_views = {thl_sec_2d, qw_sec_2d, wthl_sec_2d, wqw_sec_2d, qwthl_sec_2d, uw_sec_2d, vw_sec_2d, wtke_sec_2d, w_sec_2d};
  ekat::device_to_host({thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec}, dim1, dim2, host_2d_views);
}

void compute_brunt_shoc_length_host(Int nlev, Int nlevi, Int shcol, Real* dz_zt, Real* thv, Real* thv_zi, Real* brunt)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_2d> temp_d(4);
  std::vector<int> dim1_sizes(4, shcol);
  std::vector<int> dim2_sizes        = {nlev,  nlev,  nlevi,  nlev};
  std::vector<const Real*> ptr_array = {dz_zt, thv,   thv_zi, brunt};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    dz_zt_d (temp_d[0]),
    thv_d   (temp_d[1]),
    thv_zi_d(temp_d[2]),
    brunt_d (temp_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto dz_zt_s  = ekat::subview(dz_zt_d, i);
    const auto thv_s    = ekat::subview(thv_d, i);
    const auto thv_zi_s = ekat::subview(thv_zi_d, i);
    const auto brunt_s  = ekat::subview(brunt_d, i);

    SHF::compute_brunt_shoc_length(team, nlev, nlevi, dz_zt_s, thv_s, thv_zi_s, brunt_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {brunt_d};
  ekat::device_to_host({brunt}, shcol, nlev, inout_views);
}

void compute_l_inf_shoc_length_host(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  // Sync to device
  std::vector<view_2d> temp_d(3);
  ekat::host_to_device({zt_grid, dz_zt, tke}, shcol, nlev, temp_d);

  // inputs
  view_2d
    zt_grid_d(temp_d[0]),
    dz_zt_d  (temp_d[1]),
    tke_d    (temp_d[2]);

  // outputs
  view_1d
    l_inf_d("l_inf", shcol);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto dz_zt_s   = ekat::subview(dz_zt_d, i);
    const auto tke_s     = ekat::subview(tke_d, i);

    Scalar l_inf_s{0};

    SHF::compute_l_inf_shoc_length(team, nlev, zt_grid_s, dz_zt_s, tke_s, l_inf_s);

    l_inf_d(i) = l_inf_s;
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {l_inf_d};
  ScreamDeepCopy::copy_to_host({l_inf}, shcol, inout_views);
}

void check_length_scale_shoc_length_host(Int nlev, Int shcol, Real* host_dx, Real* host_dy, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  // Sync to device
  std::vector<view_1d> temp_1d_d(2);
  std::vector<view_2d> temp_2d_d(1);
  ScreamDeepCopy::copy_to_device({host_dx,host_dy}, shcol, temp_1d_d);
  ekat::host_to_device({shoc_mix}, shcol, nlev, temp_2d_d);

  view_1d
    host_dx_d(temp_1d_d[0]),
    host_dy_d(temp_1d_d[1]);

  view_2d
    shoc_mix_d(temp_2d_d[0]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar host_dx_s{host_dx_d(i)};
    const Scalar host_dy_s{host_dy_d(i)};
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    SHF::check_length_scale_shoc_length(team, nlev, host_dx_s, host_dy_s, shoc_mix_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {shoc_mix_d};
  ekat::device_to_host({shoc_mix}, shcol, nlev, inout_views);
}

void shoc_diag_obklen_host(Int shcol, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc, Real* thl_sfc,
                        Real* cldliq_sfc, Real* qv_sfc, Real* ustar, Real* kbfs, Real* obklen)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using view_1d    = typename SHF::view_1d<Scalar>;

  // Sync to device
  std::vector<view_1d> temp_d(7);
  std::vector<const Real*> ptr_array = {uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, thl_sfc,
                                        cldliq_sfc, qv_sfc};
  ScreamDeepCopy::copy_to_device(ptr_array, shcol, temp_d);

  // Inputs
  view_1d
    uw_sfc_d(temp_d[0]),
    vw_sfc_d(temp_d[1]),
    wthl_sfc_d(temp_d[2]),
    wqw_sfc_d(temp_d[3]),
    thl_sfc_d(temp_d[4]),
    cldliq_sfc_d(temp_d[5]),
    qv_sfc_d(temp_d[6]);

  // Outputs
  view_1d
    ustar_d("ustar", shcol),
    kbfs_d("kbfs", shcol),
    obklen_d("obklen", shcol);

  Kokkos::parallel_for("shoc_diag_obklen", shcol, KOKKOS_LAMBDA (const int& i) {
    Scalar uw_sfc_s{uw_sfc_d(i)};
    Scalar vw_sfc_s{vw_sfc_d(i)};
    Scalar wthl_sfc_s{wthl_sfc_d(i)};
    Scalar wqw_sfc_s{wqw_sfc_d(i)};
    Scalar thl_sfc_s{thl_sfc_d(i)};
    Scalar cldliq_sfc_s{cldliq_sfc_d(i)};
    Scalar qv_sfc_s{qv_sfc_d(i)};

    Scalar ustar_s{0};
    Scalar kbfs_s{0};
    Scalar obklen_s{0};

    SHF::shoc_diag_obklen(uw_sfc_s, vw_sfc_s, wthl_sfc_s, wqw_sfc_s, thl_sfc_s, cldliq_sfc_s, qv_sfc_s,
                          ustar_s, kbfs_s, obklen_s);

    ustar_d(i) = ustar_s;
    kbfs_d(i) = kbfs_s;
    obklen_d(i) = obklen_s;
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {ustar_d, kbfs_d, obklen_d};
  ScreamDeepCopy::copy_to_host({ustar, kbfs, obklen}, shcol, inout_views);
}

void shoc_pblintd_cldcheck_host(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh) {
  using SHOC    = Functions<Real, DefaultDevice>;
  using Spack   = typename SHOC::Spack;
  using Scalar  = typename SHOC::Scalar;
  using view_2d = typename SHOC::view_2d<Spack>;
  using view_1d = typename SHOC::view_1d<Scalar>;

  std::vector<Int> dim1(2, shcol);
  std::vector<Int> dim2 = {nlevi, nlev};

  std::vector<view_2d> cldcheck_2d(2);
  ekat::host_to_device({zi, cldn}, dim1, dim2, cldcheck_2d);

  view_2d
    zi_2d  (cldcheck_2d[0]),
    cldn_2d(cldcheck_2d[1]);

  std::vector<view_1d> cldcheck_1d(1);
  ScreamDeepCopy::copy_to_device({pblh}, shcol, cldcheck_1d);

  view_1d pblh_1d (cldcheck_1d[0]);

  Kokkos::parallel_for("pblintd_cldcheck", shcol, KOKKOS_LAMBDA (const int& i) {

    const int nlev_v_indx = (nlev-1)/Spack::n;
    const int nlev_p_indx = (nlev-1)%Spack::n;
    Scalar zi_s   = zi_2d  (i, nlev_v_indx)[nlev_p_indx];
    Scalar cldn_s = cldn_2d(i, nlev_v_indx)[nlev_p_indx];
    Scalar pblh_s = pblh_1d(i);

    SHOC::shoc_pblintd_cldcheck(zi_s, cldn_s, pblh_s);

    pblh_1d(i) = pblh_s;
  });

  std::vector<view_1d> inout_views = {pblh_1d};
  ScreamDeepCopy::copy_to_host({pblh}, shcol, inout_views);
}

void shoc_length_host(Int shcol, Int nlev, Int nlevi, Real* host_dx, Real* host_dy,
                   Real* zt_grid, Real* zi_grid, Real*dz_zt, Real* tke,
                   Real* thv, Real*brunt, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_1d> temp_1d_d(2);
  std::vector<view_2d> temp_2d_d(7);
  std::vector<int> dim1_sizes(7, shcol);
  std::vector<int> dim2_sizes = {nlev, nlevi, nlev, nlev,
                                 nlev, nlev, nlev};
  std::vector<const Real*> ptr_array = {zt_grid, zi_grid, dz_zt, tke,
                                        thv,     brunt,   shoc_mix};
  // Sync to device
  ScreamDeepCopy::copy_to_device({host_dx, host_dy}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);

  // inputs
  view_1d
    host_dx_d(temp_1d_d[0]),
    host_dy_d(temp_1d_d[1]);

  view_2d
    zt_grid_d(temp_2d_d[0]),
    zi_grid_d(temp_2d_d[1]),
    dz_zt_d(temp_2d_d[2]),
    tke_d(temp_2d_d[3]),
    thv_d(temp_2d_d[4]),
    brunt_d(temp_2d_d[5]),
    shoc_mix_d(temp_2d_d[6]);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 1, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    // Inputs
    const Scalar host_dx_s{host_dx_d(i)};
    const Scalar host_dy_s{host_dy_d(i)};

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto thv_s = ekat::subview(thv_d, i);
    const auto brunt_s = ekat::subview(brunt_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    // Hardcode runtime option for F90 tests.
    const Scalar length_fac = 0.5;
    SHF::shoc_length(team,nlev,nlevi,length_fac,host_dx_s,host_dy_s,
                     zt_grid_s,zi_grid_s,dz_zt_s,tke_s,
                     thv_s,workspace,brunt_s,shoc_mix_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {brunt_d,shoc_mix_d};
  ekat::device_to_host({brunt,shoc_mix}, shcol, nlev, inout_views);
}

void shoc_energy_fixer_host(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* zt_grid,
                         Real* zi_grid, Real* se_b, Real* ke_b, Real* wv_b, Real* wl_b,
                         Real* se_a, Real* ke_a, Real* wv_a, Real* wl_a, Real* wthl_sfc,
                         Real* wqw_sfc, Real* rho_zt, Real* tke, Real* pint,
                         Real* host_dse)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_1d> temp_1d_d(10);
  std::vector<const Real*> ptr_array_1d = {se_b, ke_b, wv_b, wl_b,     se_a,
                                           ke_a, wv_a, wl_a, wthl_sfc, wqw_sfc};
  std::vector<view_2d> temp_2d_d(6);
  std::vector<int> dim1_sizes(6, shcol);
  std::vector<int> dim2_sizes           = {nlev,    nlevi,   nlevi,
                                           nlev,    nlev,    nlev};
  std::vector<const Real*> ptr_array_2d = {zt_grid, zi_grid, pint,
                                           rho_zt,  tke,     host_dse};

  // Sync to device
  ScreamDeepCopy::copy_to_device(ptr_array_1d, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array_2d, dim1_sizes, dim2_sizes, temp_2d_d);

  view_1d
    se_b_d(temp_1d_d[0]),
    ke_b_d(temp_1d_d[1]),
    wv_b_d(temp_1d_d[2]),
    wl_b_d(temp_1d_d[3]),
    se_a_d(temp_1d_d[4]),
    ke_a_d(temp_1d_d[5]),
    wv_a_d(temp_1d_d[6]),
    wl_a_d(temp_1d_d[7]),
    wthl_sfc_d(temp_1d_d[8]),
    wqw_sfc_d(temp_1d_d[9]);

  view_2d
    zt_grid_d(temp_2d_d[0]),
    zi_grid_d(temp_2d_d[1]),
    pint_d(temp_2d_d[2]),
    rho_zt_d(temp_2d_d[3]),
    tke_d(temp_2d_d[4]),
    host_dse_d(temp_2d_d[5]);


  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 1, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar se_b_s{se_b_d(i)};
    const Scalar ke_b_s{ke_b_d(i)};
    const Scalar wv_b_s{wv_b_d(i)};
    const Scalar wl_b_s{wl_b_d(i)};
    const Scalar se_a_s{se_a_d(i)};
    const Scalar ke_a_s{ke_a_d(i)};
    const Scalar wv_a_s{wv_a_d(i)};
    const Scalar wl_a_s{wl_a_d(i)};
    const Scalar wthl_sfc_s{wthl_sfc_d(i)};
    const Scalar wqw_sfc_s{wqw_sfc_d(i)};

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto pint_s = ekat::subview(pint_d, i);
    const auto rho_zt_s = ekat::subview(rho_zt_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto host_dse_s = ekat::subview(host_dse_d, i);

    SHF::shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid_s,zi_grid_s,se_b_s,
                           ke_b_s,wv_b_s,wl_b_s,se_a_s,ke_a_s,wv_a_s,wl_a_s,
                           wthl_sfc_s,wqw_sfc_s,rho_zt_s,tke_s,pint_s,workspace,host_dse_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {host_dse_d};
  ekat::device_to_host({host_dse}, shcol, nlev, inout_views);
}

void compute_shoc_vapor_host(Int shcol, Int nlev, Real* qw, Real* ql, Real* qv)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 3;

  // Sync to device
  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device( {qw,  ql, qv}, shcol, nlev, temp_d);

  // Inputs/Outputs
  view_2d
    qw_d(temp_d[0]),
    ql_d(temp_d[1]),
    qv_d(temp_d[2]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto qw_s = ekat::subview(qw_d, i);
    const auto ql_s = ekat::subview(ql_d, i);
    const auto qv_s = ekat::subview(qv_d, i);

    SHF::compute_shoc_vapor(team, nlev, qw_s, ql_s, qv_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {qv_d};
  ekat::device_to_host({qv}, shcol, nlev, inout_views);
}

void update_prognostics_implicit_host(Int shcol, Int nlev, Int nlevi, Int num_tracer, Real dtime,
                                   Real* dz_zt, Real* dz_zi, Real* rho_zt, Real* zt_grid, Real* zi_grid,
                                   Real* tk, Real* tkh, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc,
                                   Real* wqw_sfc, Real* wtracer_sfc, Real* thetal, Real* qw, Real* tracer,
                                   Real* tke, Real* u_wind, Real* v_wind)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using view_3d    = typename SHF::view_3d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_2d_arrays = 13;

  std::vector<view_1d> temp_1d_d(4);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);
  std::vector<view_3d> temp_3d_d(1);

  std::vector<int> dim1_sizes(num_2d_arrays, shcol);
  std::vector<int> dim2_sizes = {nlev, nlevi, nlev,       nlev, nlevi,
                                 nlev, nlev,  num_tracer, nlev, nlev,
                                 nlev, nlev,  nlev};
  std::vector<const Real*> ptr_array = {dz_zt, dz_zi,  rho_zt,       zt_grid, zi_grid,
                                        tk,    tkh,    wtracer_sfc,  thetal,  qw,
                                        tke,   u_wind, v_wind};

  // Sync to device
  ScreamDeepCopy::copy_to_device({uw_sfc, vw_sfc, wthl_sfc, wqw_sfc}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);
  ekat::host_to_device({tracer}, shcol, nlev, num_tracer, temp_3d_d);

  view_1d
    uw_sfc_d(temp_1d_d[0]),
    vw_sfc_d(temp_1d_d[1]),
    wthl_sfc_d(temp_1d_d[2]),
    wqw_sfc_d(temp_1d_d[3]);

  view_2d
    dz_zt_d(temp_2d_d[0]),
    dz_zi_d(temp_2d_d[1]),
    rho_zt_d(temp_2d_d[2]),
    zt_grid_d(temp_2d_d[3]),
    zi_grid_d(temp_2d_d[4]),
    tk_d(temp_2d_d[5]),
    tkh_d(temp_2d_d[6]),
    wtracer_sfc_d(temp_2d_d[7]),
    thetal_d(temp_2d_d[8]),
    qw_d(temp_2d_d[9]),
    tke_d(temp_2d_d[10]),
    u_wind_d(temp_2d_d[11]),
    v_wind_d(temp_2d_d[12]);

  view_3d
    qtracers_f90_d(temp_3d_d[0]);

  // Local variables
  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // CXX version of shoc qtracers is the transpose of the fortran version.
  view_3d qtracers_cxx_d("",shcol,num_tracer,nlev_packs);

  // scalarize each view
  const auto qtracers_cxx_d_s = ekat::scalarize(qtracers_cxx_d);
  const auto qtracers_f90_d_s = ekat::scalarize(qtracers_f90_d);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_tracer), [&] (const Int& q) {
        qtracers_cxx_d_s(i,q,k) = qtracers_f90_d_s(i,k,q);
      });
    });
  });

  // Local variable workspace
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(num_tracer+3)*Spack::n;
  const int tmp_var_size = 8+n_wind_slots+n_trac_slots;
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, tmp_var_size, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar uw_sfc_s{uw_sfc_d(i)};
    const Scalar vw_sfc_s{vw_sfc_d(i)};
    const Scalar wthl_sfc_s{wthl_sfc_d(i)};
    const Scalar wqw_sfc_s{wqw_sfc_d(i)};

    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto rho_zt_s = ekat::subview(rho_zt_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto tk_s = ekat::subview(tk_d, i);
    const auto tkh_s = ekat::subview(tkh_d, i);
    const auto wtracer_sfc_s = ekat::subview(wtracer_sfc_d, i);
    const auto thetal_s = ekat::subview(thetal_d, i);
    const auto qw_s = ekat::subview(qw_d, i);
    const auto u_wind_s = ekat::subview(u_wind_d, i);
    const auto v_wind_s = ekat::subview(v_wind_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto tracer_s = Kokkos::subview(qtracers_cxx_d, i, Kokkos::ALL(), Kokkos::ALL());

    SHF::update_prognostics_implicit(team, nlev, nlevi, num_tracer, dtime,
                                     dz_zt_s, dz_zi_s, rho_zt_s, zt_grid_s,
                                     zi_grid_s, tk_s, tkh_s, uw_sfc_s, vw_sfc_s,
                                     wthl_sfc_s, wqw_sfc_s, wtracer_sfc_s,
                                     workspace,
                                     thetal_s, qw_s, tracer_s, tke_s, u_wind_s, v_wind_s);
  });

  // Transpose tracers
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_tracer), [&] (const Int& q) {
        qtracers_f90_d_s(i,k,q) = qtracers_cxx_d_s(i,q,k);
      });
    });
  });

  // Sync back to host
  std::vector<view_2d> inout_views_2d = {thetal_d, qw_d, u_wind_d, v_wind_d, tke_d};
  ekat::device_to_host({thetal, qw, u_wind, v_wind, tke}, shcol, nlev, inout_views_2d);

  std::vector<view_3d> inout_views = {qtracers_f90_d};
  ekat::device_to_host({tracer}, shcol, nlev, num_tracer, inout_views);
}

void diag_third_shoc_moments_host(Int shcol, Int nlev, Int nlevi, Real* w_sec, Real* thl_sec,
                               Real* wthl_sec, Real* isotropy, Real* brunt, Real* thetal,
                               Real* tke, Real* dz_zt, Real* dz_zi, Real* zt_grid, Real* zi_grid,
                               Real* w3)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_2d> temp_d(12);
  std::vector<int> dim1_sizes(12, shcol);
  std::vector<int> dim2_sizes = {nlev, nlevi, nlevi,  nlev,
                                 nlev,  nlev,  nlev,  nlev,
                                 nlevi, nlev,  nlevi, nlevi};
  std::vector<const Real*> ptr_array = {w_sec, thl_sec, wthl_sec, isotropy,
                                        brunt, thetal,  tke,      dz_zt,
                                        dz_zi, zt_grid, zi_grid,  w3};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    wsec_d(temp_d[0]),
    thl_sec_d(temp_d[1]),
    wthl_sec_d(temp_d[2]),
    isotropy_d(temp_d[3]),
    brunt_d(temp_d[4]),
    thetal_d(temp_d[5]),
    tke_d(temp_d[6]),
    dz_zt_d(temp_d[7]),
    dz_zi_d(temp_d[8]),
    zt_grid_d(temp_d[9]),
    zi_grid_d(temp_d[10]),
    w3_d(temp_d[11]);

  // Local variables
  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 4, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const auto wsec_s = ekat::subview(wsec_d, i);
    const auto thl_sec_s = ekat::subview(thl_sec_d, i);
    const auto wthl_sec_s = ekat::subview(wthl_sec_d, i);
    const auto isotropy_s = ekat::subview(isotropy_d, i);
    const auto brunt_s = ekat::subview(brunt_d, i);
    const auto thetal_s = ekat::subview(thetal_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto w3_s = ekat::subview(w3_d, i);

    // Hardcode for F90 testing
    const Real c_diag_3rd_mom = 7.0;
    SHF::diag_third_shoc_moments(team, nlev, nlevi, c_diag_3rd_mom, wsec_s, thl_sec_s,
                                 wthl_sec_s, isotropy_s, brunt_s, thetal_s, tke_s,
                                 dz_zt_s, dz_zi_s, zt_grid_s, zi_grid_s,
                                 workspace,
                                 w3_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {w3_d};
  ekat::device_to_host({w3}, shcol, nlevi, inout_views);
}

void adv_sgs_tke_host(Int nlev, Int shcol, Real dtime, Real* shoc_mix, Real* wthv_sec,
                   Real* sterm_zt, Real* tk, Real* tke, Real* a_diss)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 6;

  std::vector<view_2d> temp_d(num_arrays);
  std::vector<const Real*> ptr_array  = {shoc_mix, wthv_sec, sterm_zt, tk,    tke,   a_diss};

  // Sync to device
  ekat::host_to_device(ptr_array, shcol, nlev, temp_d);

  view_2d
    //input
    shoc_mix_d (temp_d[0]),
    wthv_sec_d (temp_d[1]),
    sterm_zt_d (temp_d[2]),
    tk_d       (temp_d[3]),
    //output
    tke_d      (temp_d[4]), //inout
    a_diss_d   (temp_d[5]); //out

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();

      // inputs
      const auto shoc_mix_s = ekat::subview(shoc_mix_d ,i);
      const auto wthv_sec_s = ekat::subview(wthv_sec_d ,i);
      const auto sterm_zt_s = ekat::subview(sterm_zt_d ,i);
      const auto tk_s       = ekat::subview(tk_d ,i);
      const auto tke_s      = ekat::subview(tke_d ,i);
      const auto a_diss_s   = ekat::subview(a_diss_d ,i);

      SHF::adv_sgs_tke(team, nlev, dtime, shoc_mix_s, wthv_sec_s, sterm_zt_s, tk_s, tke_s, a_diss_s);
    });

  // Sync back to host
  std::vector<view_2d> inout_views = {tke_d, a_diss_d};
  ekat::device_to_host({tke, a_diss}, shcol, nlev, inout_views);
}

void shoc_assumed_pdf_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* w_field,
                        Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* w_sec, Real* wqw_sec,
                        Real* qwthl_sec, Real* w3, Real* pres, Real* zt_grid, Real* zi_grid,
                        Real* shoc_cldfrac, Real* shoc_ql, Real* wqls, Real* wthv_sec, Real* shoc_ql2)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 18;

  std::vector<view_2d> temp_d(num_arrays);
  std::vector<int> dim1_sizes(num_arrays, shcol);
  std::vector<int> dim2_sizes = {nlev,  nlev,  nlevi, nlevi, nlevi, nlev,
                                 nlevi, nlevi, nlevi, nlev,  nlev,  nlev,
                                 nlevi, nlev,  nlev,  nlev,  nlev,  nlev};
  std::vector<const Real*> ptr_array = {thetal,  qw,           thl_sec, qw_sec,  wthl_sec, w_sec,
                                        wqw_sec, qwthl_sec,    w3,      w_field, pres,     zt_grid,
                                        zi_grid, shoc_cldfrac, shoc_ql, wqls,    wthv_sec, shoc_ql2};
  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  // Inputs/Outputs
  view_2d
    thetal_d(temp_d[0]),
    qw_d(temp_d[1]),
    thl_sec_d(temp_d[2]),
    qw_sec_d(temp_d[3]),
    wthl_sec_d(temp_d[4]),
    w_sec_d(temp_d[5]),
    wqw_sec_d(temp_d[6]),
    qwthl_sec_d(temp_d[7]),
    w3_d(temp_d[8]),
    w_field_d(temp_d[9]),
    pres_d(temp_d[10]),
    zt_grid_d(temp_d[11]),
    zi_grid_d(temp_d[12]),
    shoc_cldfrac_d(temp_d[13]),
    shoc_ql_d(temp_d[14]),
    wqls_d(temp_d[15]),
    wthv_sec_d(temp_d[16]),
    shoc_ql2_d(temp_d[17]);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlev_packs, 6, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const auto thetal_s = ekat::subview(thetal_d, i);
    const auto qw_s = ekat::subview(qw_d, i);
    const auto thl_sec_s = ekat::subview(thl_sec_d, i);
    const auto qw_sec_s = ekat::subview(qw_sec_d, i);
    const auto wthl_sec_s = ekat::subview(wthl_sec_d, i);
    const auto w_sec_s = ekat::subview(w_sec_d, i);
    const auto wqw_sec_s = ekat::subview(wqw_sec_d, i);
    const auto qwthl_sec_s = ekat::subview(qwthl_sec_d, i);
    const auto w3_s = ekat::subview(w3_d, i);
    const auto w_field_s = ekat::subview(w_field_d, i);
    const auto pres_s = ekat::subview(pres_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto shoc_cldfrac_s = ekat::subview(shoc_cldfrac_d, i);
    const auto shoc_ql_s = ekat::subview(shoc_ql_d, i);
    const auto wqls_s = ekat::subview(wqls_d, i);
    const auto wthv_sec_s = ekat::subview(wthv_sec_d, i);
    const auto shoc_ql2_s = ekat::subview(shoc_ql2_d, i);

    SHF::shoc_assumed_pdf(team, nlev, nlevi, thetal_s, qw_s, w_field_s, thl_sec_s, qw_sec_s, wthl_sec_s, w_sec_s,
                          wqw_sec_s, qwthl_sec_s, w3_s, pres_s, zt_grid_s, zi_grid_s,
                          workspace,
                          shoc_cldfrac_s, shoc_ql_s, wqls_s, wthv_sec_s, shoc_ql2_s);
  });

  // Sync back to host
  std::vector<view_2d> out_views = {shoc_cldfrac_d, shoc_ql_d, wqls_d, wthv_sec_d, shoc_ql2_d};
  ekat::device_to_host({shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2}, shcol, nlev, out_views);
}
void compute_shr_prod_host(Int nlevi, Int nlev, Int shcol, Real* dz_zi, Real* u_wind, Real* v_wind, Real* sterm)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 4;

  std::vector<view_2d> temp_d(num_arrays);
  std::vector<int> dim1_sizes(num_arrays, shcol);
  std::vector<int> dim2_sizes = {nlevi,   nlev,   nlev, nlevi};
  std::vector<const Real*> ptr_array  = {dz_zi, u_wind, v_wind, sterm};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  view_2d
    //input
    dz_zi_d (temp_d[0]),
    u_wind_d(temp_d[1]),
    v_wind_d(temp_d[2]),
    //output
    sterm_d (temp_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();

      // inputs
      const auto dz_zi_s  = ekat::subview(dz_zi_d ,i);
      const auto u_wind_s = ekat::subview(u_wind_d ,i);
      const auto v_wind_s = ekat::subview(v_wind_d ,i);
      //output
      const auto sterm_s  = ekat::subview(sterm_d ,i);

      SHF::compute_shr_prod(team, nlevi, nlev, dz_zi_s, u_wind_s, v_wind_s, sterm_s);
    });

  // Sync back to host
  std::vector<view_2d> inout_views = {sterm_d};
  ekat::device_to_host({sterm}, shcol, nlevi, inout_views);
}

void compute_tmpi_host(Int nlevi, Int shcol, Real dtime, Real *rho_zi, Real *dz_zi, Real *tmpi)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 3;

  // Sync to device
  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device({rho_zi,  dz_zi, tmpi}, shcol, nlevi, temp_d);

  // Inputs/Outputs
  view_2d
    rho_zi_d(temp_d[0]),
    dz_zi_d(temp_d[1]),
    tmpi_d(temp_d[2]);

  const Int nk_pack = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto rho_zi_s = ekat::subview(rho_zi_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto tmpi_s = ekat::subview(tmpi_d, i);

    SHF::compute_tmpi(team, nlevi, dtime, rho_zi_s, dz_zi_s, tmpi_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tmpi_d};
  ekat::device_to_host({tmpi}, shcol, nlevi, inout_views);
}

void integ_column_stability_host(Int nlev, Int shcol, Real *dz_zt,
                              Real *pres, Real* brunt, Real *brunt_int)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 3;

  // Sync to device
  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device({dz_zt, pres, brunt}, shcol, nlev, temp_d);

  // Inputs
  view_2d
    dz_zt_d(temp_d[0]),
    pres_d (temp_d[1]),
    brunt_d(temp_d[2]);

  //Output
  view_1d brunt_int_d("brunt_int", shcol);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();

      // create subviews of the views
      const auto dz_zt_s = ekat::subview(dz_zt_d, i);
      const auto pres_s  = ekat::subview(pres_d, i);
      const auto brunt_s = ekat::subview(brunt_d, i);

      //declare output as a scalar
      Scalar brunt_int_s{0};

      SHF::integ_column_stability(team, nlev, dz_zt_s, pres_s, brunt_s, brunt_int_s);

      brunt_int_d(i) = brunt_int_s;
    });

  // Sync back to host
  std::vector<view_1d> inout_views = {brunt_int_d};
  ScreamDeepCopy::copy_to_host({brunt_int}, shcol, inout_views);
}

void isotropic_ts_host(Int nlev, Int shcol, Real* brunt_int, Real* tke,
                    Real* a_diss, Real* brunt, Real* isotropy)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  std::vector<view_1d> temp_1d(1); // for 1d array

  static constexpr Int num_arrays = 4;
  std::vector<view_2d> temp_2d(num_arrays); //for 2d arrays
  std::vector<const Real*> ptr_array = {tke, a_diss, brunt, isotropy};

  // Sync to device
  ScreamDeepCopy::copy_to_device({brunt_int}, shcol, temp_1d);
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d);

  //inputs
  view_1d brunt_int_d(temp_1d[0]);

  view_2d
    tke_d     (temp_2d[0]),
    a_diss_d  (temp_2d[1]),
    brunt_d   (temp_2d[2]),
    //output
    isotropy_d(temp_2d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();

      // inputs
      const Scalar brunt_int_s{brunt_int_d(i)};
      const auto tke_s     = ekat::subview(tke_d,      i);// create subviews of the views
      const auto a_diss_s  = ekat::subview(a_diss_d,   i);
      const auto brunt_s   = ekat::subview(brunt_d,    i);

      //outputs
      const auto isotropy_s = ekat::subview(isotropy_d, i); //output

      // Hard code these runtime options for F90
      const Real lambda_low = 0.001;
      const Real lambda_high   = 0.08;
      const Real lambda_slope  = 2.65;
      const Real lambda_thresh = 0.02;
      SHF::isotropic_ts(team, nlev, lambda_low, lambda_high, lambda_slope, lambda_thresh,
		      brunt_int_s, tke_s, a_diss_s, brunt_s, isotropy_s);
    });

  // Sync back to host
  std::vector<view_2d> inout_views = {isotropy_d};
  ekat::device_to_host({isotropy}, shcol, nlev, inout_views);

}

void dp_inverse_host(Int nlev, Int shcol, Real *rho_zt, Real *dz_zt, Real *rdp_zt)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 3;

  // Sync to device
  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device({rho_zt,  dz_zt, rdp_zt}, shcol, nlev, temp_d);

  // Inputs/Outputs
  view_2d
    rho_zt_d(temp_d[0]),
    dz_zt_d(temp_d[1]),
    rdp_zt_d(temp_d[2]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto rho_zt_s = ekat::subview(rho_zt_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto rdp_zt_s = ekat::subview(rdp_zt_d, i);

    SHF::dp_inverse(team, nlev, rho_zt_s, dz_zt_s, rdp_zt_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {rdp_zt_d};
  ekat::device_to_host({rdp_zt}, shcol, nlev, inout_views);
}

int shoc_init_host(Int nlev, Real *pref_mid, Int nbot_shoc, Int ntop_shoc)
{
  using SHF  = Functions<Real, DefaultDevice>;
  using Spack       = typename SHF::Spack;
  using view_1d     = typename SHF::view_1d<Spack>;

  // Sync to device
  std::vector<view_1d> temp_d(1);
  ekat::host_to_device({pref_mid}, nlev, temp_d);
  view_1d pref_mid_d(temp_d[0]);

  return SHF::shoc_init(nbot_shoc,ntop_shoc,pref_mid_d);
}

Int shoc_main_host(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Int npbl, Real* host_dx, Real* host_dy, Real* thv, Real* zt_grid,
                Real* zi_grid, Real* pres, Real* presi, Real* pdel, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc,
                Real* wtracer_sfc, Int num_qtracers, Real* w_field, Real* inv_exner, Real* phis, Real* host_dse, Real* tke,
                Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk,
                Real* shoc_ql, Real* shoc_cldfrac, Real* pblh, Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec,
                Real* qw_sec, Real* qwthl_sec, Real* wthl_sec, Real* wqw_sec, Real* wtke_sec, Real* uw_sec, Real* vw_sec,
                Real* w3, Real* wqls_sec, Real* brunt, Real* shoc_ql2)
{

  using SHF  = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using view_3d    = typename SHF::view_3d<Spack>;
  using ExeSpace   = typename SHF::KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  // Initialize Kokkos views, sync to device
  static constexpr Int num_1d_arrays = 7;
  static constexpr Int num_2d_arrays = 35;
  static constexpr Int num_3d_arrays = 1;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);
  std::vector<view_3d> temp_3d_d(num_3d_arrays);

  std::vector<int> dim1_2d_sizes = {shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol};
  std::vector<int> dim2_2d_sizes = {nlev,  nlevi, nlev,         nlevi, nlev,
                                    nlev,  nlev,  num_qtracers, nlev,  nlev,
                                    nlev,  nlev,  nlev,         nlev,  nlev,
                                    nlev,  nlev,  nlev,         nlev,  nlev,
                                    nlev,  nlev,  nlev,         nlevi, nlevi,
                                    nlevi, nlevi, nlevi,        nlevi, nlevi,
                                    nlevi, nlevi, nlev,         nlev,  nlev};

  std::vector<const Real*> ptr_array_1d = {host_dx, host_dy, wthl_sfc, wqw_sfc,
                                           uw_sfc,  vw_sfc,  phis};
  std::vector<const Real*> ptr_array_2d = {zt_grid,   zi_grid,  pres,          presi,        pdel,
                                           thv,       w_field,  wtracer_sfc,   inv_exner,    host_dse,
                                           tke,       thetal,   qw,            u_wind,       v_wind,
                                           wthv_sec,  tk,       shoc_cldfrac,  shoc_ql,      shoc_ql2,
                                           tkh,       shoc_mix, w_sec,         thl_sec,      qw_sec,
                                           qwthl_sec, wthl_sec, wqw_sec,       wtke_sec,     uw_sec,
                                           vw_sec,    w3,       wqls_sec,      brunt,        isotropy};

  ScreamDeepCopy::copy_to_device(ptr_array_1d, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array_2d, dim1_2d_sizes, dim2_2d_sizes, temp_2d_d);
  ekat::host_to_device({qtracers}, shcol, nlev, num_qtracers, temp_3d_d);

  Int index_counter = 0;
  view_1d
    host_dx_d (temp_1d_d[index_counter++]),
    host_dy_d (temp_1d_d[index_counter++]),
    wthl_sfc_d(temp_1d_d[index_counter++]),
    wqw_sfc_d (temp_1d_d[index_counter++]),
    uw_sfc_d  (temp_1d_d[index_counter++]),
    vw_sfc_d  (temp_1d_d[index_counter++]),
    phis_d    (temp_1d_d[index_counter++]),
    pblh_d    ("pblh",shcol),
    ustar_d   ("ustar",shcol),
    obklen_d  ("obklen",shcol);

  index_counter = 0;
  view_2d
    zt_grid_d     (temp_2d_d[index_counter++]),
    zi_grid_d     (temp_2d_d[index_counter++]),
    pres_d        (temp_2d_d[index_counter++]),
    presi_d       (temp_2d_d[index_counter++]),
    pdel_d        (temp_2d_d[index_counter++]),
    thv_d         (temp_2d_d[index_counter++]),
    w_field_d     (temp_2d_d[index_counter++]),
    wtracer_sfc_d (temp_2d_d[index_counter++]),
    inv_exner_d   (temp_2d_d[index_counter++]),
    host_dse_d    (temp_2d_d[index_counter++]),
    tke_d         (temp_2d_d[index_counter++]),
    thetal_d      (temp_2d_d[index_counter++]),
    qw_d          (temp_2d_d[index_counter++]),
    u_wind_d      (temp_2d_d[index_counter++]),
    v_wind_d      (temp_2d_d[index_counter++]),
    wthv_sec_d    (temp_2d_d[index_counter++]),
    tk_d          (temp_2d_d[index_counter++]),
    shoc_cldfrac_d(temp_2d_d[index_counter++]),
    shoc_ql_d     (temp_2d_d[index_counter++]),
    shoc_ql2_d    (temp_2d_d[index_counter++]),
    tkh_d         (temp_2d_d[index_counter++]),
    shoc_mix_d    (temp_2d_d[index_counter++]),
    w_sec_d       (temp_2d_d[index_counter++]),
    thl_sec_d     (temp_2d_d[index_counter++]),
    qw_sec_d      (temp_2d_d[index_counter++]),
    qwthl_sec_d   (temp_2d_d[index_counter++]),
    wthl_sec_d    (temp_2d_d[index_counter++]),
    wqw_sec_d     (temp_2d_d[index_counter++]),
    wtke_sec_d    (temp_2d_d[index_counter++]),
    uw_sec_d      (temp_2d_d[index_counter++]),
    vw_sec_d      (temp_2d_d[index_counter++]),
    w3_d          (temp_2d_d[index_counter++]),
    wqls_sec_d    (temp_2d_d[index_counter++]),
    brunt_d       (temp_2d_d[index_counter++]),
    isotropy_d    (temp_2d_d[index_counter++]);

  view_3d
    qtracers_f90_d(temp_3d_d[0]);

  // shoc_main treats u/v_wind as 1 array and
  // CXX version of shoc qtracers is the transpose of the fortran version..
  const auto nlev_packs = ekat::npack<Spack>(nlev);
  view_3d horiz_wind_d("horiz_wind",shcol,2,nlev_packs);
  view_3d qtracers_cxx_d("qtracers",shcol,num_qtracers,nlev_packs);

  // scalarize each view
  const auto u_wind_d_s = ekat::scalarize(u_wind_d);
  const auto v_wind_d_s = ekat::scalarize(v_wind_d);
  const auto horiz_wind_d_s = ekat::scalarize(horiz_wind_d);
  const auto qtracers_cxx_d_s = ekat::scalarize(qtracers_cxx_d);
  const auto qtracers_f90_d_s = ekat::scalarize(qtracers_f90_d);

  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
      horiz_wind_d_s(i,0,k) = u_wind_d_s(i,k);
      horiz_wind_d_s(i,1,k) = v_wind_d_s(i,k);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
        qtracers_cxx_d_s(i,q,k) = qtracers_f90_d_s(i,k,q);
      });
    });
  });

  // Pack our data into structs and ship it off to shoc_main.
  SHF::SHOCInput shoc_input{host_dx_d,  host_dy_d,     zt_grid_d,   zi_grid_d,
                             pres_d,    presi_d,       pdel_d,      thv_d,
                             w_field_d, wthl_sfc_d,    wqw_sfc_d,   uw_sfc_d,
                             vw_sfc_d,  wtracer_sfc_d, inv_exner_d, phis_d};
  SHF::SHOCInputOutput shoc_input_output{host_dse_d,   tke_d,      thetal_d,       qw_d,
                                         horiz_wind_d, wthv_sec_d, qtracers_cxx_d,
                                         tk_d,         shoc_cldfrac_d, shoc_ql_d};
  SHF::SHOCOutput shoc_output{pblh_d, ustar_d, obklen_d, shoc_ql2_d, tkh_d};
  SHF::SHOCHistoryOutput shoc_history_output{shoc_mix_d,  w_sec_d,    thl_sec_d, qw_sec_d,
                                             qwthl_sec_d, wthl_sec_d, wqw_sec_d, wtke_sec_d,
                                             uw_sec_d,    vw_sec_d,   w3_d,      wqls_sec_d,
                                             brunt_d,     isotropy_d};
  SHF::SHOCRuntime shoc_runtime_options{0.001,0.04,2.65,0.02,1.0,1.0,1.0,1.0,0.5,7.0,0.1,0.1};

  const auto nlevi_packs = ekat::npack<Spack>(nlevi);

#ifdef SCREAM_SHOC_SMALL_KERNELS
  view_1d
    se_b   ("se_b", shcol),
    ke_b   ("ke_b", shcol),
    wv_b   ("wv_b", shcol),
    wl_b   ("wl_b", shcol),
    se_a   ("se_a", shcol),
    ke_a   ("ke_a", shcol),
    wv_a   ("wv_a", shcol),
    wl_a   ("wl_a", shcol),
    kbfs   ("kbfs", shcol),
    ustar2 ("ustar2", shcol),
    wstar  ("wstar", shcol);

  view_2d
    rho_zt  ("rho_zt",  shcol, nlevi_packs),
    shoc_qv ("shoc_qv", shcol, nlevi_packs),
    tabs    ("shoc_tabs", shcol, nlev_packs),
    dz_zt   ("dz_zt",   shcol, nlevi_packs),
    dz_zi   ("dz_zi",   shcol, nlevi_packs);

  SHF::SHOCTemporaries shoc_temporaries{
    se_b, ke_b, wv_b, wl_b, se_a, ke_a, wv_a, wl_a, kbfs, ustar2, wstar,
    rho_zt, shoc_qv, tabs, dz_zt, dz_zi};
#endif

  // Create local workspace
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
  ekat::WorkspaceManager<Spack, SHF::KT::Device> workspace_mgr(nlevi_packs, 14+(n_wind_slots+n_trac_slots), policy);

  const auto elapsed_microsec = SHF::shoc_main(shcol, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
                                               workspace_mgr, shoc_runtime_options,
                                               shoc_input, shoc_input_output, shoc_output, shoc_history_output
#ifdef SCREAM_SHOC_SMALL_KERNELS
                                               , shoc_temporaries
#endif
                                               );

  // Copy wind back into separate views and
  // Transpose tracers
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
      u_wind_d_s(i,k) = horiz_wind_d_s(i,0,k);
      v_wind_d_s(i,k) = horiz_wind_d_s(i,1,k);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
        qtracers_f90_d_s(i,k,q) = qtracers_cxx_d_s(i,q,k);
      });
    });
  });

  // Sync back to host
  // 1d
  std::vector<view_1d> out_views_1d = {pblh_d};
  ScreamDeepCopy::copy_to_host({pblh}, shcol, out_views_1d);

  // 2d
  std::vector<int> dim1_2d_out = {shcol, shcol, shcol, shcol, shcol,
                                  shcol, shcol, shcol, shcol,
                                  shcol, shcol, shcol, shcol, shcol,
                                  shcol, shcol, shcol, shcol, shcol,
                                  shcol, shcol, shcol, shcol, shcol,
                                  shcol};
  std::vector<int> dim2_2d_out = {nlev,  nlev,  nlev,  nlev,  nlev,
                                  nlev,  nlev,  nlev,  nlev,
                                  nlev,  nlev,  nlev,  nlev,  nlevi,
                                  nlevi, nlevi, nlevi, nlevi, nlevi,
                                  nlevi, nlevi, nlevi, nlev,  nlev,
                                  nlev};
  std::vector<Real*> ptr_array_2d_out = {host_dse, tke,       thetal,   qw,       u_wind,
                                         v_wind,   wthv_sec,  tk,       shoc_cldfrac,
                                         shoc_ql,  shoc_ql2,  shoc_mix, w_sec,    thl_sec,
                                         qw_sec,   qwthl_sec, wthl_sec, wqw_sec,  wtke_sec,
                                         uw_sec,   vw_sec,    w3,       wqls_sec, brunt,
                                         isotropy};
  std::vector<view_2d> out_views_2d = {host_dse_d, tke_d,       thetal_d,   qw_d,       u_wind_d,
                                       v_wind_d,   wthv_sec_d,  tk_d,       shoc_cldfrac_d,
                                       shoc_ql_d,  shoc_ql2_d,  shoc_mix_d, w_sec_d,    thl_sec_d,
                                       qw_sec_d,   qwthl_sec_d, wthl_sec_d, wqw_sec_d,  wtke_sec_d,
                                       uw_sec_d,   vw_sec_d,    w3_d,       wqls_sec_d, brunt_d,
                                       isotropy_d};
  ekat::device_to_host(ptr_array_2d_out, dim1_2d_out, dim2_2d_out, out_views_2d);

  // 3d
  std::vector<view_3d> out_views_3d = {qtracers_f90_d};
  ekat::device_to_host({qtracers}, shcol, nlev, num_qtracers, out_views_3d);

  return elapsed_microsec;
}

void pblintd_height_host(Int shcol, Int nlev, Int npbl, Real* z, Real* u, Real* v, Real* ustar, Real* thv, Real* thv_ref, Real* pblh, Real* rino, bool* check)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using Scalar     = typename SHOC::Scalar;
  using view_1d    = typename SHOC::view_1d<Scalar>;
  using bview_1d   = typename SHOC::view_1d<bool>;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using ExeSpace   = typename SHOC::KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  std::vector<view_2d> views_2d(5);
  ekat::host_to_device({z, u, v, thv, rino}, shcol, nlev, views_2d);

  view_2d z_2d   (views_2d[0]),
          u_2d   (views_2d[1]),
          v_2d   (views_2d[2]),
          thv_2d (views_2d[3]),
          rino_2d(views_2d[4]);

  std::vector<view_1d> views_1d(3);
  ScreamDeepCopy::copy_to_device({ustar, thv_ref, pblh}, shcol, views_1d);
  view_1d ustar_1d   (views_1d[0]),
          thv_ref_1d (views_1d[1]),
          pblh_1d    (views_1d[2]);

  std::vector<bview_1d> views_bool_1d(1);
  ScreamDeepCopy::copy_to_device({check}, shcol, views_bool_1d);
  bview_1d check_1d (views_bool_1d[0]);

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto z_1d    = ekat::subview(z_2d, i);
    const auto u_1d    = ekat::subview(u_2d, i);
    const auto v_1d    = ekat::subview(v_2d, i);
    const auto thv_1d  = ekat::subview(thv_2d, i);
    const auto rino_1d = ekat::subview(rino_2d, i);

    Scalar& ustar_s   = ustar_1d(i);
    Scalar& thv_ref_s = thv_ref_1d(i);
    Scalar& pblh_s    = pblh_1d(i);
    bool& check_s     = check_1d(i);

    SHOC::pblintd_height(team, nlev, npbl, z_1d, u_1d, v_1d, ustar_s, thv_1d, thv_ref_s, pblh_s, rino_1d, check_s);
 });

  std::vector<view_1d> out_1d_views = {pblh_1d};
  ScreamDeepCopy::copy_to_host({pblh}, shcol, out_1d_views);

  std::vector<view_2d> out_2d_views = {rino_2d};
  ekat::device_to_host({rino}, shcol, nlev, out_2d_views);

  std::vector<bview_1d> out_bool_1d_views = {check_1d};
  ScreamDeepCopy::copy_to_host({check}, shcol, out_bool_1d_views);
}

void vd_shoc_decomp_and_solve_host(Int shcol, Int nlev, Int nlevi, Int num_rhs, Real dtime, Real* kv_term, Real* tmpi, Real* rdp_zt, Real* flux, Real* var)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar         = typename SHF::Scalar;
  using Spack          = typename SHF::Spack;
  using view_1d        = typename SHF::view_1d<Scalar>;
  using view_2d        = typename SHF::view_2d<Spack>;
  using view_2d_scalar = typename SHF::view_2d<Scalar>;
  using view_3d        = typename SHF::view_3d<Spack>;
  using KT             = typename SHF::KT;
  using ExeSpace       = typename KT::ExeSpace;
  using MemberType     = typename SHF::MemberType;

  static constexpr Int num_1d_arrays = 1;
  static constexpr Int num_2d_arrays = 3;
  static constexpr Int num_3d_arrays = 1;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);
  std::vector<view_3d> temp_3d_d(num_3d_arrays);

  std::vector<int> dim1_sizes(num_2d_arrays, shcol);
  std::vector<int> dim2_sizes = {nlevi, nlevi, nlev};
  std::vector<const Real*> ptr_array = {kv_term, tmpi, rdp_zt};

  // Sync to device
  ScreamDeepCopy::copy_to_device({flux}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);
  ekat::host_to_device({var}, shcol, nlev, num_rhs, temp_3d_d);

  view_1d
    flux_d(temp_1d_d[0]);

  view_2d
    kv_term_d(temp_2d_d[0]),
    tmpi_d(temp_2d_d[1]),
    rdp_zt_d(temp_2d_d[2]);

  view_2d_scalar
    du_d("du", shcol, nlev),
    dl_d("dl", shcol, nlev),
    d_d ("d",  shcol, nlev);

  view_3d
    var_d(temp_3d_d[0]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar flux_s{flux_d(i)};
    const auto kv_term_s = ekat::subview(kv_term_d, i);
    const auto tmpi_s    = ekat::subview(tmpi_d, i);
    const auto rdp_zt_s  = ekat::subview(rdp_zt_d, i);
    const auto du_s      = ekat::subview(du_d, i);
    const auto dl_s      = ekat::subview(dl_d, i);
    const auto d_s       = ekat::subview(d_d, i);
    const auto var_s = Kokkos::subview(var_d, i, Kokkos::ALL(), Kokkos::ALL());

    SHF::vd_shoc_decomp(team, nlev, kv_term_s, tmpi_s, rdp_zt_s, dtime, flux_s, du_s, dl_s, d_s);
    team.team_barrier();
    SHF::vd_shoc_solve(team, du_s, dl_s, d_s, var_s);
  });

  // Sync back to host
  std::vector<view_3d> inout_views = {var_d};
  ekat::device_to_host({var}, shcol, nlev, num_rhs, inout_views);
}

void shoc_grid_host(Int shcol, Int nlev, Int nlevi, Real* zt_grid, Real* zi_grid, Real* pdel, Real* dz_zt, Real* dz_zi, Real* rho_zt)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_2d_arrays = 6;
  std::vector<view_2d> temp_2d_d(num_2d_arrays);
  std::vector<int> dim1_sizes(num_2d_arrays, shcol);
  std::vector<int> dim2_sizes = {nlev, nlevi, nlev,
                                 nlev, nlevi, nlev};
  std::vector<const Real*> ptr_array = {zt_grid, zi_grid, pdel,
                                        dz_zt,   dz_zi,   rho_zt};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);

  view_2d
    zt_grid_d(temp_2d_d[0]),
    zi_grid_d(temp_2d_d[1]),
    pdel_d(temp_2d_d[2]),
    dz_zt_d(temp_2d_d[3]),
    dz_zi_d(temp_2d_d[4]),
    rho_zt_d(temp_2d_d[5]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto pdel_s = ekat::subview(pdel_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto rho_zt_s = ekat::subview(rho_zt_d, i);

    SHF::shoc_grid(team, nlev, nlevi, zt_grid_s, zi_grid_s, pdel_s, dz_zt_s, dz_zi_s, rho_zt_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {dz_zt_d, dz_zi_d, rho_zt_d};
  ekat::device_to_host<Int>({dz_zt, dz_zi, rho_zt}, {shcol, shcol, shcol}, {nlev, nlevi, nlev}, inout_views);
}

void eddy_diffusivities_host(Int nlev, Int shcol, Real* pblh, Real* zt_grid, Real* tabs, Real* shoc_mix, Real* sterm_zt,
                          Real* isotropy, Real* tke, Real* tkh, Real* tk)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_1d_arrays = 1;
  static constexpr Int num_2d_arrays = 8;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);

  std::vector<const Real*> ptr_array = {zt_grid,  tabs, shoc_mix, sterm_zt,
                                        isotropy, tke,  tkh,      tk};

  // Sync to device
  ScreamDeepCopy::copy_to_device({pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d);

  view_1d pblh_d(temp_1d_d[0]);

  view_2d
    zt_grid_d(temp_2d_d[0]),
    tabs_d(temp_2d_d[1]),
    shoc_mix_d(temp_2d_d[2]),
    sterm_zt_d(temp_2d_d[3]),
    isotropy_d(temp_2d_d[4]),
    tke_d(temp_2d_d[5]),
    tkh_d(temp_2d_d[6]),
    tk_d(temp_2d_d[7]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar pblh_s{pblh_d(i)};

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto tabs_s = ekat::subview(tabs_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);
    const auto sterm_zt_s = ekat::subview(sterm_zt_d, i);
    const auto isotropy_s = ekat::subview(isotropy_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto tkh_s = ekat::subview(tkh_d, i);
    const auto tk_s = ekat::subview(tk_d, i);

    // Hardcode runtime options for F90 testing
    const Real Ckh = 0.1;
    const Real Ckm = 0.1;
    SHF::eddy_diffusivities(team, nlev, Ckh, Ckm, pblh_s, zt_grid_s, tabs_s, shoc_mix_s, sterm_zt_s, isotropy_s, tke_s, tkh_s, tk_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tkh_d, tk_d};
  ekat::device_to_host({tkh, tk}, shcol, nlev, inout_views);
}

void pblintd_surf_temp_host(Int shcol, Int nlev, Int nlevi, Real* z, Real* ustar, Real* obklen, Real* kbfs, Real* thv, Real* tlv, Real* pblh, bool* check, Real* rino)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using Scalar     = typename SHOC::Scalar;
  using view_1d      = typename SHOC::view_1d<Scalar>;
  using view_bool_1d = typename SHOC::view_1d<bool>;
  using view_2d      = typename SHOC::view_2d<Spack>;

  std::vector<view_2d> views_2d(3);
  ekat::host_to_device({z, thv, rino}, shcol, nlev, views_2d);
  view_2d z_2d   (views_2d[0]),
          thv_2d (views_2d[1]),
          rino_2d(views_2d[2]);

  std::vector<view_1d> views_1d(5);
  ScreamDeepCopy::copy_to_device({ustar, obklen, kbfs, pblh, tlv}, shcol, views_1d);
  view_1d ustar_1d (views_1d[0]),
          obklen_1d(views_1d[1]),
          kbfs_1d  (views_1d[2]),
          pblh_1d  (views_1d[3]),
          tlv_1d   (views_1d[4]);

  std::vector<view_bool_1d> bviews_1d(1);
  ScreamDeepCopy::copy_to_device({check}, shcol, bviews_1d);
  view_bool_1d check_1d (bviews_1d[0]);

  Int npbl = nlev;

  Kokkos::parallel_for("pblintd_surf_temp", shcol, KOKKOS_LAMBDA (const int& i) {

    const auto z_1d    = ekat::subview(z_2d, i);
    const auto thv_1d  = ekat::subview(thv_2d, i);
    const auto rino_1d = ekat::subview(rino_2d, i);

    const Scalar& ustar_s  = ustar_1d(i);
    const Scalar& obklen_s = obklen_1d(i);
    const Scalar& kbfs_s   = kbfs_1d(i);
    Scalar& pblh_s         = pblh_1d(i);

    Scalar& tlv_s = tlv_1d(i);
    bool& check_s = check_1d(i);

    SHOC::pblintd_surf_temp(nlev, nlevi, npbl, z_1d, ustar_s, obklen_s, kbfs_s,
               thv_1d, tlv_s, pblh_s, check_s, rino_1d);
  });

  std::vector<view_1d> out_1d_views = {pblh_1d, tlv_1d};
  ScreamDeepCopy::copy_to_host({pblh, tlv}, shcol, out_1d_views);

  std::vector<view_2d> out_2d_views = {rino_2d};
  ekat::device_to_host({rino}, shcol, nlev, out_2d_views);

  std::vector<view_bool_1d> out_bool_1d_views = {check_1d};
  ScreamDeepCopy::copy_to_host({check}, shcol, out_bool_1d_views);
}

void pblintd_check_pblh_host(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* ustar, bool* check, Real* pblh)
{
  using SHOC         = Functions<Real, DefaultDevice>;
  using Spack        = typename SHOC::Spack;
  using Scalar       = typename SHOC::Scalar;
  using view_bool_1d = typename SHOC::view_1d<bool>;
  using view_1d      = typename SHOC::view_1d<Scalar>;
  using view_2d      = typename SHOC::view_2d<Spack>;

  std::vector<view_2d> views_2d(1);
  ekat::host_to_device({z}, shcol, nlev, views_2d);
  view_2d z_2d (views_2d[0]);

  std::vector<view_1d> views_1d(2);
  ScreamDeepCopy::copy_to_device({ustar, pblh}, shcol, views_1d);
  view_1d ustar_1d (views_1d[0]),
          pblh_1d  (views_1d[1]);

  std::vector<view_bool_1d> bool_views_1d(1);
  ScreamDeepCopy::copy_to_device({check}, shcol, bool_views_1d);
  view_bool_1d check_1d(bool_views_1d[0]);

  Kokkos::parallel_for("pblintd_check_pblh", shcol, KOKKOS_LAMBDA (const int& i) {

    const auto z_1d  = ekat::subview(z_2d, i);
    Scalar& ustar_s  = ustar_1d(i);
    Scalar& pblh_s   = pblh_1d(i);
    bool check_s     = (bool)(check_1d(i));

    SHOC::pblintd_check_pblh(nlevi, npbl, z_1d, ustar_s, check_s, pblh_s);
 });

  std::vector<view_1d> out_1d_views = {pblh_1d};
  ScreamDeepCopy::copy_to_host({pblh}, shcol, out_1d_views);
}

void pblintd_host(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* zi, Real* thl, Real* ql, Real* q, Real* u, Real* v, Real* ustar, Real* obklen, Real* kbfs, Real* cldn, Real* pblh)
{
    using SHF = Functions<Real, DefaultDevice>;

    using Scalar     = typename SHF::Scalar;
    using Spack      = typename SHF::Spack;
    using view_1d    = typename SHF::view_1d<Scalar>;
    using view_2d    = typename SHF::view_2d<Spack>;
    using KT         = typename SHF::KT;
    using ExeSpace   = typename KT::ExeSpace;
    using MemberType = typename SHF::MemberType;

    static constexpr Int num_1d_arrays = 3;
    static constexpr Int num_2d_arrays = 8;

    std::vector<view_1d> temp_1d_d(num_1d_arrays);
    std::vector<view_2d> temp_2d_d(num_2d_arrays);

    std::vector<int> dim1_sizes(num_2d_arrays, shcol);
    std::vector<int> dim2_sizes = {nlev, nlevi, nlev, nlev,
                                   nlev, nlev,  nlev, nlev};
    std::vector<const Real*> ptr_array = {z, zi, thl, ql,
                                          q, u,  v,   cldn};

    // Sync to device
    ScreamDeepCopy::copy_to_device({ustar, obklen, kbfs}, shcol, temp_1d_d);
    ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);

    view_1d
      ustar_d(temp_1d_d[0]),
      obklen_d(temp_1d_d[1]),
      kbfs_d(temp_1d_d[2]),
      pblh_d("pblh", shcol);

    view_2d
      z_d(temp_2d_d[0]),
      zi_d(temp_2d_d[1]),
      thl_d(temp_2d_d[2]),
      ql_d(temp_2d_d[3]),
      q_d(temp_2d_d[4]),
      u_d(temp_2d_d[5]),
      v_d(temp_2d_d[6]),
      cldn_d(temp_2d_d[7]);

    const Int nlev_pack = ekat::npack<Spack>(nlev);
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_pack);

    // Local variable workspace
    ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlev_pack, 2, policy);

    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();

      auto workspace = workspace_mgr.get_workspace(team);

      const Scalar ustar_s{ustar_d(i)};
      const Scalar obklen_s{obklen_d(i)};
      const Scalar kbfs_s{kbfs_d(i)};
      Scalar pblh_s{0};

      const auto z_s = ekat::subview(z_d, i);
      const auto zi_s = ekat::subview(zi_d, i);
      const auto thl_s = ekat::subview(thl_d, i);
      const auto ql_s = ekat::subview(ql_d, i);
      const auto q_s = ekat::subview(q_d, i);
      const auto u_s = ekat::subview(u_d, i);
      const auto v_s = ekat::subview(v_d, i);
      const auto cldn_s = ekat::subview(cldn_d, i);

      SHF::pblintd(team,nlev,nlevi,npbl,z_s,zi_s,thl_s,ql_s,q_s,u_s,v_s,ustar_s,
                           obklen_s,kbfs_s,cldn_s,workspace,pblh_s);

      pblh_d(i) = pblh_s;
    });

    // Sync back to host
    std::vector<view_1d> out_views = {pblh_d};
    ScreamDeepCopy::copy_to_host({pblh}, shcol, out_views);
}

void shoc_tke_host(Int shcol, Int nlev, Int nlevi, Real dtime, Real* wthv_sec, Real* shoc_mix, Real* dz_zi, Real* dz_zt, Real* pres,
                Real* tabs, Real* u_wind, Real* v_wind, Real* brunt, Real* zt_grid, Real* zi_grid, Real* pblh, Real* tke, Real* tk,
                Real* tkh, Real* isotropy)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using view_1d    = typename SHF::view_1d<Scalar>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_1d_arrays = 1;
  static constexpr Int num_2d_arrays = 15;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);

  std::vector<int> dim1_sizes(num_2d_arrays, shcol);
  std::vector<int> dim2_sizes = {nlev, nlev, nlev,  nlev, nlevi, nlev, nlev, nlev,
                                 nlev, nlev, nlevi, nlev, nlev,  nlev, nlev};
  std::vector<const Real*> ptr_array = {wthv_sec, shoc_mix, u_wind,  v_wind, dz_zi, dz_zt, pres, tabs,
                                        brunt,    zt_grid,  zi_grid, tke,    tk,    tkh,   isotropy};

  // Sync to device
  ScreamDeepCopy::copy_to_device({pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d);

  view_1d pblh_d(temp_1d_d[0]);

  view_2d
    wthv_sec_d(temp_2d_d[0]),
    shoc_mix_d(temp_2d_d[1]),
    u_wind_d(temp_2d_d[2]),
    v_wind_d(temp_2d_d[3]),
    dz_zi_d(temp_2d_d[4]),
    dz_zt_d(temp_2d_d[5]),
    pres_d(temp_2d_d[6]),
    tabs_d(temp_2d_d[7]),
    brunt_d(temp_2d_d[8]),
    zt_grid_d(temp_2d_d[9]),
    zi_grid_d(temp_2d_d[10]),
    tke_d(temp_2d_d[11]),
    tk_d(temp_2d_d[12]),
    tkh_d(temp_2d_d[13]),
    isotropy_d(temp_2d_d[14]);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 3, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar pblh_s{pblh_d(i)};

    const auto wthv_sec_s = ekat::subview(wthv_sec_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);
    const auto u_wind_s = ekat::subview(u_wind_d, i);
    const auto v_wind_s = ekat::subview(v_wind_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto pres_s = ekat::subview(pres_d, i);
    const auto tabs_s = ekat::subview(tabs_d, i);
    const auto brunt_s = ekat::subview(brunt_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto tk_s = ekat::subview(tk_d, i);
    const auto tkh_s = ekat::subview(tkh_d, i);
    const auto isotropy_s = ekat::subview(isotropy_d, i);

    // Hardcode for F90 testing
    const Real lambda_low    = 0.001;
    const Real lambda_high   = 0.08;
    const Real lambda_slope  = 2.65;
    const Real lambda_thresh = 0.02;
    const Real Ckh           = 0.1;
    const Real Ckm           = 0.1;

    SHF::shoc_tke(team,nlev,nlevi,dtime,lambda_low,lambda_high,lambda_slope,lambda_thresh,
                  Ckh, Ckm,
		  wthv_sec_s,shoc_mix_s,dz_zi_s,dz_zt_s,pres_s,
                  tabs_s,u_wind_s,v_wind_s,brunt_s,zt_grid_s,zi_grid_s,pblh_s,
                  workspace,
                  tke_s,tk_s,tkh_s,isotropy_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tke_d, tk_d, tkh_d, isotropy_d};
  ekat::device_to_host({tke, tk, tkh, isotropy}, shcol, nlev, inout_views);
}

void compute_shoc_temperature_host(Int shcol, Int nlev, Real *thetal, Real *ql, Real *inv_exner, Real* tabs)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 4;

  // Sync to device
  std::vector<view_2d> temp_d(num_arrays);
  ekat::host_to_device({thetal, ql, inv_exner, tabs}, shcol, nlev, temp_d);

  // Inputs/Outputs
  view_2d
    thetal_d(temp_d[0]),
    ql_d(temp_d[1]),
    inv_exner_d(temp_d[2]),
    tabs_d(temp_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto thetal_s = ekat::subview(thetal_d, i);
    const auto ql_s = ekat::subview(ql_d, i);
    const auto inv_exner_s = ekat::subview(inv_exner_d, i);
    const auto tabs_s = ekat::subview(tabs_d, i);

    SHF::compute_shoc_temperature(team, nlev, thetal_s, ql_s, inv_exner_s, tabs_s);
  });

  // Sync back to host
  std::vector<view_2d> out_views = {tabs_d};
  ekat::device_to_host({tabs}, shcol, nlev, out_views);
}

void shoc_assumed_pdf_tilde_to_real_host(Real w_first, Real sqrtw2, Real* w1)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_w1(*w1);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack sqrtw2_(sqrtw2), w1_(local_w1), w_first_(w_first);
    SHF::shoc_assumed_pdf_tilde_to_real(w_first_, sqrtw2_, w1_);
    t_d(0) = w1_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *w1 = t_h(0);
}

void shoc_assumed_pdf_vv_parameters_host(Real w_first, Real w_sec, Real w3var, Real w_tol_sqd, Real* skew_w, Real* w1_1, Real* w1_2, Real* w2_1, Real* w2_2, Real* a)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 6);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack w3var_(w3var), w_first_(w_first), w_sec_(w_sec), a_, skew_w_, w1_1_, w1_2_, w2_1_, w2_2_;
    SHF::shoc_assumed_pdf_vv_parameters(w_first_, w_sec_, w3var_, w_tol_sqd, skew_w_, w1_1_, w1_2_, w2_1_, w2_2_, a_);
    t_d(0) = a_[0];
    t_d(1) = skew_w_[0];
    t_d(2) = w1_1_[0];
    t_d(3) = w1_2_[0];
    t_d(4) = w2_1_[0];
    t_d(5) = w2_2_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *a = t_h(0);
  *skew_w = t_h(1);
  *w1_1 = t_h(2);
  *w1_2 = t_h(3);
  *w2_1 = t_h(4);
  *w2_2 = t_h(5);
}

void shoc_assumed_pdf_thl_parameters_host(Real wthlsec, Real sqrtw2, Real sqrtthl, Real thlsec, Real thl_first, Real w1_1, Real w1_2, Real skew_w, Real a, Real thl_tol, Real w_thresh, Real* thl1_1, Real* thl1_2, Real* thl2_1, Real* thl2_2, Real* sqrtthl2_1, Real* sqrtthl2_2)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 6);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), skew_w_(skew_w), sqrtthl_(sqrtthl), sqrtw2_(sqrtw2), thl_first_(thl_first), thlsec_(thlsec), w1_1_(w1_1), w1_2_(w1_2), wthlsec_(wthlsec), sqrtthl2_1_, sqrtthl2_2_, thl1_1_, thl1_2_, thl2_1_, thl2_2_;
    SHF::shoc_assumed_pdf_thl_parameters(wthlsec_, sqrtw2_, sqrtthl_, thlsec_, thl_first_, w1_1_, w1_2_, skew_w_, a_, thl_tol, w_thresh, thl1_1_, thl1_2_, thl2_1_, thl2_2_, sqrtthl2_1_, sqrtthl2_2_);
    t_d(0) = sqrtthl2_1_[0];
    t_d(1) = sqrtthl2_2_[0];
    t_d(2) = thl1_1_[0];
    t_d(3) = thl1_2_[0];
    t_d(4) = thl2_1_[0];
    t_d(5) = thl2_2_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *sqrtthl2_1 = t_h(0);
  *sqrtthl2_2 = t_h(1);
  *thl1_1 = t_h(2);
  *thl1_2 = t_h(3);
  *thl2_1 = t_h(4);
  *thl2_2 = t_h(5);
}

void shoc_assumed_pdf_qw_parameters_host(Real wqwsec, Real sqrtw2, Real skew_w, Real sqrtqt, Real qwsec, Real w1_2, Real w1_1, Real qw_first, Real a, Real rt_tol, Real w_thresh, Real* qw1_1, Real* qw1_2, Real* qw2_1, Real* qw2_2, Real* sqrtqw2_1, Real* sqrtqw2_2)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 6);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), qw_first_(qw_first), qwsec_(qwsec), skew_w_(skew_w), sqrtqt_(sqrtqt), sqrtw2_(sqrtw2), w1_1_(w1_1), w1_2_(w1_2), wqwsec_(wqwsec), qw1_1_, qw1_2_, qw2_1_, qw2_2_, sqrtqw2_1_, sqrtqw2_2_;
    SHF::shoc_assumed_pdf_qw_parameters(wqwsec_, sqrtw2_, skew_w_, sqrtqt_, qwsec_, w1_2_, w1_1_, qw_first_, a_, rt_tol, w_thresh, qw1_1_, qw1_2_, qw2_1_, qw2_2_, sqrtqw2_1_, sqrtqw2_2_);
    t_d(0) = qw1_1_[0];
    t_d(1) = qw1_2_[0];
    t_d(2) = qw2_1_[0];
    t_d(3) = qw2_2_[0];
    t_d(4) = sqrtqw2_1_[0];
    t_d(5) = sqrtqw2_2_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *qw1_1 = t_h(0);
  *qw1_2 = t_h(1);
  *qw2_1 = t_h(2);
  *qw2_2 = t_h(3);
  *sqrtqw2_1 = t_h(4);
  *sqrtqw2_2 = t_h(5);
}

void shoc_assumed_pdf_inplume_correlations_host(Real sqrtqw2_1, Real sqrtthl2_1, Real a, Real sqrtqw2_2, Real sqrtthl2_2, Real qwthlsec, Real qw1_1, Real qw_first, Real thl1_1, Real thl_first, Real qw1_2, Real thl1_2, Real* r_qwthl_1)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), qw1_1_(qw1_1), qw1_2_(qw1_2), qw_first_(qw_first), qwthlsec_(qwthlsec), sqrtqw2_1_(sqrtqw2_1), sqrtqw2_2_(sqrtqw2_2), sqrtthl2_1_(sqrtthl2_1), sqrtthl2_2_(sqrtthl2_2), thl1_1_(thl1_1), thl1_2_(thl1_2), thl_first_(thl_first), r_qwthl_1_;
    SHF::shoc_assumed_pdf_inplume_correlations(sqrtqw2_1_, sqrtthl2_1_, a_, sqrtqw2_2_, sqrtthl2_2_, qwthlsec_, qw1_1_, qw_first_, thl1_1_, thl_first_, qw1_2_, thl1_2_, r_qwthl_1_);
    t_d(0) = r_qwthl_1_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *r_qwthl_1 = t_h(0);
}

void shoc_assumed_pdf_compute_temperature_host(Real thl1, Real pval, Real* tl1)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack pval_(pval), thl1_(thl1), tl1_;
    SHF::shoc_assumed_pdf_compute_temperature(thl1_, pval_, tl1_);
    t_d(0) = tl1_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *tl1 = t_h(0);
}

void shoc_assumed_pdf_compute_qs_host(Real tl1_1, Real tl1_2, Real pval, Real* qs1, Real* beta1, Real* qs2, Real* beta2)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using Smask   = typename SHF::Smask;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 4);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack pval_(pval), tl1_1_(tl1_1), tl1_2_(tl1_2), beta1_, beta2_, qs1_, qs2_;
    Smask active_entries(true);
    SHF::shoc_assumed_pdf_compute_qs(tl1_1_, tl1_2_, pval_, active_entries, qs1_, beta1_, qs2_, beta2_);
    t_d(0) = beta1_[0];
    t_d(1) = beta2_[0];
    t_d(2) = qs1_[0];
    t_d(3) = qs2_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *beta1 = t_h(0);
  *beta2 = t_h(1);
  *qs1 = t_h(2);
  *qs2 = t_h(3);
}

void shoc_assumed_pdf_compute_s_host(Real qw1, Real qs1, Real beta, Real pval, Real thl2, Real qw2, Real sqrtthl2, Real sqrtqw2, Real r_qwthl, Real* s, Real* std_s, Real* qn, Real* c)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 4);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack beta_(beta), pval_(pval), qs1_(qs1), qw1_(qw1), qw2_(qw2), r_qwthl_(r_qwthl), sqrtqw2_(sqrtqw2), sqrtthl2_(sqrtthl2), thl2_(thl2), c_, qn_, s_, std_s_;
    SHF::shoc_assumed_pdf_compute_s(qw1_, qs1_, beta_, pval_, thl2_, qw2_, sqrtthl2_, sqrtqw2_, r_qwthl_, s_, std_s_, qn_, c_);
    t_d(0) = c_[0];
    t_d(1) = qn_[0];
    t_d(2) = s_[0];
    t_d(3) = std_s_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *c = t_h(0);
  *qn = t_h(1);
  *s = t_h(2);
  *std_s = t_h(3);
}

void shoc_assumed_pdf_compute_sgs_liquid_host(Real a, Real ql1, Real ql2, Real* shoc_ql)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), ql1_(ql1), ql2_(ql2), shoc_ql_;
    SHF::shoc_assumed_pdf_compute_sgs_liquid(a_, ql1_, ql2_, shoc_ql_);
    t_d(0) = shoc_ql_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *shoc_ql = t_h(0);
}

void shoc_assumed_pdf_compute_cloud_liquid_variance_host(Real a, Real s1, Real ql1, Real c1, Real std_s1, Real s2, Real ql2, Real c2, Real std_s2, Real shoc_ql, Real* shoc_ql2)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), c1_(c1), c2_(c2), ql1_(ql1), ql2_(ql2), s1_(s1), s2_(s2), shoc_ql_(shoc_ql), std_s1_(std_s1), std_s2_(std_s2), shoc_ql2_;
    SHF::shoc_assumed_pdf_compute_cloud_liquid_variance(a_, s1_, ql1_, c1_, std_s1_, s2_, ql2_, c2_, std_s2_, shoc_ql_, shoc_ql2_);
    t_d(0) = shoc_ql2_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *shoc_ql2 = t_h(0);
}

void shoc_assumed_pdf_compute_liquid_water_flux_host(Real a, Real w1_1, Real w_first, Real ql1, Real w1_2, Real ql2, Real* wqls)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack   = typename SHF::Spack;
  using view_1d = typename SHF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack a_(a), ql1_(ql1), ql2_(ql2), w1_1_(w1_1), w1_2_(w1_2), w_first_(w_first), wqls_;
    SHF::shoc_assumed_pdf_compute_liquid_water_flux(a_, w1_1_, w_first_, ql1_, w1_2_, ql2_, wqls_);
    t_d(0) = wqls_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *wqls = t_h(0);
}

void shoc_assumed_pdf_compute_buoyancy_flux_host(Real wthlsec, Real wqwsec, Real pval, Real wqls, Real* wthv_sec)
{
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;
  using view_1d = typename PF::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack pval_(pval), wqls_(wqls), wqwsec_(wqwsec), wthlsec_(wthlsec), wthv_sec_;
    PF::shoc_assumed_pdf_compute_buoyancy_flux(wthlsec_, wqwsec_, pval_, wqls_, wthv_sec_);
    t_d(0) = wthv_sec_[0];
  });
  Kokkos::deep_copy(t_h, t_d);
  *wthv_sec = t_h(0);
}

} // namespace shoc
} // namespace scream
