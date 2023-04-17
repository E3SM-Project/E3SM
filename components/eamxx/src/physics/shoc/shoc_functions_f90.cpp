#include "shoc_functions_f90.hpp"

#include "shoc_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include "share/util/scream_deep_copy.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to SHOC fortran calls. The stubs below will link to fortran definitions in shoc_iso_c.f90
//

extern "C" {

// Special shoc_init function for shoc_main_bfb test
void shoc_init_for_main_bfb_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                              Real zvir, Real latvap, Real latice, Real karman,
                              Real* pref_mid, int nbot_shoc, int ntop_shoc);
void shoc_use_cxx_c(bool use_cxx);


void shoc_grid_c(int shcol, int nlev, int nlevi, Real *zt_grid, Real *zi_grid,
                 Real *pdel, Real *dz_zt, Real *dzi_zi, Real *rho_zt);

void shoc_diag_obklen_c(Int shcol, Real *uw_sfc, Real *vw_sfc, Real *wthl_sfc,
                        Real *wqw_sfc, Real *thl_sfc, Real *cldliq_sfc,
                        Real *qv_sfc, Real *ustar, Real *kbfs, Real *obklen);

void update_host_dse_c(Int shcol, Int nlev, Real *thlm, Real *shoc_ql,
                       Real *inv_exner, Real *zt_grid, Real *phis, Real *host_dse);

void shoc_energy_fixer_c(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv,
                         Real *zt_grid, Real *zi_grid, Real *se_b, Real *ke_b,
                         Real *wv_b, Real *wl_b, Real *se_a, Real *ke_a,
                         Real *wv_a, Real *wl_a, Real *wthl_sfc, Real *wqw_sfc,
                         Real *rho_zt, Real *tke, Real *pint,
                         Real *host_dse);

void shoc_energy_integrals_c(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int);

void shoc_energy_total_fixer_c(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv,
                               Real *zt_grid, Real *zi_grid,
                               Real *se_b, Real *ke_b, Real *wv_b, Real *wl_b,
                               Real *se_a, Real *ke_a, Real *wv_a, Real *wl_a,
                               Real *wthl_sfc, Real *wqw_sfc, Real *rho_zt,
                               Real *te_a, Real *te_b);

void shoc_energy_threshold_fixer_c(Int shcol, Int nlev, Int nlevi,
                             Real *pint, Real *tke, Real *te_a, Real *te_b,
                             Real *se_dis, Int *shoctop);

void shoc_energy_dse_fixer_c(Int shcol, Int nlev,
                             Real *se_dis, Int *shoctop,
                             Real *host_dse);

void calc_shoc_varorcovar_c(Int shcol, Int nlev, Int nlevi,  Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
                            Real *invar1, Real *invar2, Real *varorcovar);

void compute_tmpi_c(Int nlevi, Int shcol, Real dtime, Real *rho_zi,
                    Real *dz_zi, Real *tmpi);

void dp_inverse_c(Int nlev, Int shcol, Real *rho_zt, Real *dz_zt, Real *rdp_zt);

void sfc_fluxes_c(Int shcol, Int num_tracer, Real dtime, Real *rho_zi_sfc,
                  Real *rdp_zt_sfc, Real *wthl_sfc, Real *wqw_sfc, Real *wtracer_sfc,
                  Real *wtke_sfc, Real *thetal, Real *qw, Real *tke, Real *tracer);

void impli_srf_stress_term_c(Int shcol, Real *rho_zi_sfc, Real *uw_sfc,
                             Real *vw_sfc, Real *u_wind_sfc, Real *v_wind_sfc,
                             Real *ksrf);

void tke_srf_flux_term_c(Int shcol, Real *uw_sfc, Real *vw_sfc,
                         Real *wtke_sfc);

void check_tke_c(Int shcol, Int nlev, Real *tke);

void shoc_tke_c(Int shcol, Int nlev, Int nlevi, Real dtime, Real *wthv_sec,
                Real *shoc_mix, Real *dz_zi, Real *dz_zt, Real *pres,
                Real *u_wind, Real *v_wind, Real *brunt, Real *obklen,
                Real *zt_grid, Real *zi_grid, Real *pblh, Real *tke,
                Real *tk, Real *tkh, Real *isotropy);

void integ_column_stability_c(Int nlev, Int shcol, Real *dz_zt, Real *pres,
                              Real *brunt, Real *brunt_int);

void compute_shr_prod_c(Int nlevi, Int nlev, Int shcol, Real *dz_zi,
                        Real *u_wind, Real *v_wind, Real *sterm);

void isotropic_ts_c(Int nlev, Int shcol, Real *brunt_int, Real *tke,
                    Real *a_diss, Real *brunt, Real *isotropy);

void adv_sgs_tke_c(Int nlev, Int shcol, Real dtime, Real *shoc_mix,
                   Real *wthv_sec, Real *sterm_zt, Real *tk,
                   Real *tke, Real *a_diss);

void eddy_diffusivities_c(Int nlev, Int shcol, Real *obklen, Real *pblh,
                          Real *zt_grid, Real *shoc_mix, Real *sterm_zt,
                          Real *isotropy, Real *tke, Real *tkh, Real *tk);

void calc_shoc_vertflux_c(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
                          Real *dz_zi, Real *invar, Real *vertflux);

void shoc_length_c(Int shcol, Int nlev, Int nlevi, Real *host_dx,
                   Real *host_dy, Real *zt_grid,
                   Real *zi_grid, Real *dz_zt,  Real *tke,
                   Real *thv, Real *brunt, Real *shoc_mix);

void compute_brunt_shoc_length_c(Int nlev, Int nlevi, Int shcol ,Real *dz_zt,
                                 Real *thv, Real *thv_zi, Real *brunt);

void compute_l_inf_shoc_length_c(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);

void compute_shoc_mix_shoc_length_c(Int nlev, Int shcol, Real *tke, Real* brunt,
                                    Real *zt_grid, Real *l_inf, Real *shoc_mix);

void check_length_scale_shoc_length_c(Int nlev, Int shcol, Real *host_dx,
                                    Real *host_dy, Real *shoc_mix);

void clipping_diag_third_shoc_moments_c(Int nlevi, Int shcol, Real *w_sec_zi,
                                    Real *w3);

void fterms_input_for_diag_third_shoc_moment_c(Real dz_zi, Real dz_zt, Real dz_zt_kc,
                                    Real isotropy_zi, Real brunt_zi, Real thetal_zi,
                                    Real *thedz, Real *thedz2, Real *iso,
                                    Real *isosqrd, Real *buoy_sgs2, Real *bet2);
void f0_to_f5_diag_third_shoc_moment_c(Real thedz, Real thedz2, Real bet2, Real iso,
                                    Real isosqrd, Real wthl_sec, Real wthl_sec_kc,
                                    Real wthl_sec_kb, Real thl_sec_kc,
                                    Real thl_sec_kb, Real w_sec, Real w_sec_kc, Real w_sec_zi,
                                    Real tke, Real tke_kc, Real *f0, Real *f1,
                                    Real *f2, Real *f3, Real *f4, Real *f5);

void omega_terms_diag_third_shoc_moment_c(Real buoy_sgs2, Real f3, Real f4,
                                    Real *omega0, Real *omega1, Real *omega2);

void x_y_terms_diag_third_shoc_moment_c(Real buoy_sgs2, Real f0, Real f1, Real f2,
                                        Real *x0, Real *y0, Real *x1, Real *y1);

void aa_terms_diag_third_shoc_moment_c(Real omega0, Real omega1, Real omega2,
                                       Real x0, Real x1, Real y0, Real y1,
                                       Real *aa0, Real *aa1);

void w3_diag_third_shoc_moment_c(Real aa0, Real aa1, Real x0,
                                 Real x1, Real f5, Real *w3);
void shoc_diag_second_moments_srf_c(Int shcol, Real* wthl_sfc, Real* uw_sfc, Real* vw_sfc,
                                    Real* ustar2, Real* wstar);

void diag_third_shoc_moments_c(Int shoc, Int nlev, Int nlevi, Real *w_sec,
                               Real *thl_sec,
                               Real *wthl_sec, Real *isotropy, Real *brunt,
                               Real *thetal, Real *tke,
                               Real *dz_zt, Real *dz_zi, Real *zt_grid,
                               Real *zi_grid, Real *w3);

void compute_diag_third_shoc_moment_c(Int shcol, Int nlev, Int nlevi, Real *w_sec,
                                      Real *thl_sec, Real *wthl_sec, Real *tke,
                                      Real *dz_zt, Real *dz_zi, Real *isotropy_zi,
                                      Real *brunt_zi, Real *w_sec_zi, Real *thetal_zi,
                                      Real *w3);

void linear_interp_c(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh);

void shoc_assumed_pdf_c(Int shcol, Int nlev, Int nlevi, Real *thetal, Real *qw,
                        Real *w_first, Real *thl_sec, Real *qw_sec, Real *wthl_sec,
                        Real *w_sec, Real *wqw_sec, Real *qwthl_sec, Real *w3,
                        Real *pres, Real *zt_grid, Real *zi_grid,
                        Real *shoc_cldfrac, Real *shoc_ql, Real *wqls,
                        Real *wthv_sec, Real *shoc_ql2);

void shoc_assumed_pdf_tilde_to_real_c(Real w_first, Real sqrtw2, Real* w1);

void shoc_assumed_pdf_vv_parameters_c(Real w_first, Real w_sec, Real w3var,
                                      Real *Skew_w, Real *w1_1, Real *w1_2,
                                      Real *w2_1, Real *w2_2, Real *a);

void shoc_assumed_pdf_thl_parameters_c(Real wthlsec, Real sqrtw2, Real sqrtthl,
                                       Real thlsec, Real thl_first, Real w1_1,
                                       Real w1_2, Real Skew_w, Real a, bool dothetal_skew,
                                       Real *thl1_1, Real *thl1_2, Real *thl2_1,
                                       Real *thl2_2, Real *sqrtthl2_1,
                                       Real *sqrtthl2_2);

void shoc_assumed_pdf_qw_parameters_c(Real wqwsec, Real sqrtw2, Real Skew_w,
                                      Real sqrtqt, Real qw_sec, Real w1_1,
                                      Real w1_2, Real qw_first, Real a,
                                      Real *qw1_1, Real *qw1_2, Real *qw2_1,
                                      Real *qw2_2, Real *sqrtqw2_1,
                                      Real *sqrtqw2_2);

void shoc_assumed_pdf_inplume_correlations_c(Real sqrtqw2_1, Real sqrtthl2_1,
                                             Real a, Real sqrtqw2_2, Real sqrtthl2_2,
                                             Real qwthlsec, Real qw1_1, Real qw_first,
                                             Real thl1_1, Real thl_first, Real qw1_2,
                                             Real thl1_2, Real *r_qwthl_1);

void shoc_assumed_pdf_compute_temperature_c(Real thl1, Real basepres,
                                            Real pval, Real *Tl1);

void shoc_assumed_pdf_compute_qs_c(Real Tl1_1, Real Tl1_2, Real pval,
                                   Real *qs1, Real *beta1, Real *qs2, Real *beta2);

void shoc_assumed_pdf_compute_s_c(Real qw1, Real qs1, Real beta, Real pval, Real thl2,
                                  Real qw2,Real sqrtthl2, Real sqrtqw2, Real r_qwthl,
                                  Real *s, Real *std_s, Real *qn, Real *C);

void shoc_assumed_pdf_compute_sgs_liquid_c(Real a, Real ql1, Real ql2, Real *shoc_ql);

void shoc_assumed_pdf_compute_cloud_liquid_variance_c(Real a, Real s1, Real ql1,
                                                      Real C1, Real std_s1, Real s2, Real ql2, Real C2,
                                                      Real std_s2, Real shoc_ql, Real *shoc_ql2);

void shoc_assumed_pdf_compute_liquid_water_flux_c(Real a, Real w1_1, Real w_first,
                                                  Real ql1, Real w1_2, Real ql2, Real *wqls);

void shoc_assumed_pdf_compute_buoyancy_flux_c(Real wthlsec, Real epsterm, Real wqwsec,
                                              Real pval, Real wqls, Real *wthv_sec);

void shoc_diag_second_moments_ubycond_c(Int shcol, Real* thl_sec, Real* qw_sec,
                                        Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec,
                                        Real* uw_sec, Real* vw_sec, Real* wtke_sec);

void shoc_pblintd_init_pot_c(Int shcol, Int nlev, Real* thl, Real* ql, Real* q, Real* thv);

void diag_second_moments_lbycond_c(Int shcol, Real *wthl_sfc, Real *wqw_sfc, Real *uw_sfc,
                                   Real *vw_sfc, Real *ustar2, Real *wstar, Real *wthl_sec,
                                   Real *wqw_sec, Real *uw_sec, Real *vw_sec, Real *wtke_sec,
                                   Real *thl_sec, Real *qw_sec, Real *qwthl_sec);

void diag_second_moments_c(Int shcol, Int nlev, Int nlevi, Real *thetal, Real *qw,
                           Real *u_wind, Real *v_wind, Real *tke, Real *isotropy,
                           Real *tkh, Real *tk, Real *dz_zi, Real *zt_grid, Real *zi_grid,
                           Real *shoc_mix, Real *thl_sec, Real *qw_sec, Real *wthl_sec,
                           Real *wqw_sec, Real *qwthl_sec, Real *uw_sec, Real *vw_sec,
                           Real *wtke_sec, Real *w_sec);

void diag_second_shoc_moments_c(Int shcol, Int nlev, Int nlevi, Real *thetal,
                                Real *qw, Real *u_wind, Real *v_wind, Real *tke,
                                Real *isotropy, Real *tkh, Real *tk, Real *dz_zi,
                                Real *zt_grid, Real *zi_grid, Real *shoc_mix,
                                Real *wthl_sfc, Real *wqw_sfc, Real *uw_sfc,
                                Real *vw_sfc, Real *thl_sec, Real *qw_sec,
                                Real *wthl_sec, Real *wqw_sec, Real *qwthl_sec,
                                Real *uw_sec, Real *vw_sec, Real *wtke_sec, Real *w_sec);

void shoc_pblintd_cldcheck_c(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh);

void compute_shoc_vapor_c(Int shcol, Int nlev, Real* qw, Real* ql, Real* qv);

void update_prognostics_implicit_c(Int shcol, Int nlev, Int nlevi, Int num_tracer, Real dtime,
                                   Real* dz_zt, Real* dz_zi, Real* rho_zt, Real* zt_grid, Real* zi_grid,
                                   Real* tk, Real* tkh, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc,
                                   Real* wqw_sfc, Real* wtracer_sfc, Real* thetal, Real* qw, Real* tracer,
                                   Real* tke, Real* u_wind, Real* v_wind);

void shoc_main_c(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* host_dx, Real* host_dy,
                 Real* thv, Real* zt_grid, Real* zi_grid, Real* pres, Real* presi, Real* pdel,
                 Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* wtracer_sfc,
                 Int num_qtracers, Real* w_field, Real* inv_exner, Real* phis, Real* host_dse, Real* tke,
                 Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* qtracers, Real* wthv_sec,
                 Real* tkh, Real* tk, Real* shoc_ql, Real* shoc_cldfrac, Real* pblh, Real* shoc_mix,
                 Real* isotropy, Real* w_sec, Real* thl_sec, Real* qw_sec, Real* qwthl_sec, Real* wthl_sec,
                 Real* wqw_sec, Real* wtke_sec, Real* uw_sec, Real* vw_sec, Real* w3, Real* wqls_sec,
                 Real* brunt, Real* shoc_ql2, Real* elapsed_s);

void pblintd_height_c(Int shcol, Int nlev, Int npbl_in, Real* z, Real* u, Real* v, Real* ustar, Real* thv, Real* thv_ref, Real* pblh, Real* rino, bool* check);

void vd_shoc_decomp_c(Int shcol, Int nlev, Int nlevi, Real* kv_term, Real* tmpi, Real* rdp_zt, Real dtime,
                      Real* flux, Real* du, Real* dl, Real* d);

void vd_shoc_solve_c(Int shcol, Int nlev, Real* du, Real* dl, Real* d, Real* var);
void pblintd_surf_temp_c(Int shcol, Int nlev, Int nlevi, Real* z, Real* ustar, Real* obklen, Real* kbfs, Real* thv, Real* tlv, Real* pblh, bool* check, Real* rino);
void pblintd_check_pblh_c(Int shcol, Int nlev, Int nlevi, Real* z, Real* ustar, bool* check, Real* pblh);
void pblintd_c(Int shcol, Int nlev, Int nlevi, Int npbl_in, Real* z, Real* zi, Real* thl, Real* ql, Real* q, Real* u, Real* v, Real* ustar, Real* obklen, Real* kbfs, Real* cldn, Real* pblh);
} // extern "C" : end _c decls

namespace scream {
namespace shoc {

//
// Glue functions to call fortran from from C++ with the Data struct
//
// In all C++ -> Fortran bridge functions you should see shoc_init(nlev, true).
// We are provisionally following P3 here in case SHOC uses global data. The
// 'true' argument is to set shoc to use its fortran implementations instead of
// calling back to C++. We want this behavior since it doesn't make much sense
// for C++ to bridge over to fortran only to have fortran bridge back to C++.
// Anyone who wants the C++ implementation should call it directly. We need
// need to be aware of data layout since f90 is different from cxx. All these
// functions will expect incoming data to be C layout. They will transpose to f90
// before calling fortran and then back to C before returning.
//

void shoc_grid(ShocGridData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt, d.dz_zi, d.rho_zt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_diag_obklen(ShocDiagObklenData& d)
{
  shoc_init(1, true); // single level function
  shoc_diag_obklen_c(d.shcol, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc, d.thl_sfc, d.cldliq_sfc, d.qv_sfc, d.ustar, d.kbfs, d.obklen);
}

void update_host_dse(UpdateHostDseData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  update_host_dse_c(d.shcol, d.nlev, d.thlm, d.shoc_ql, d.inv_exner, d.zt_grid, d.phis, d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_fixer(ShocEnergyFixerData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_fixer_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.zt_grid, d.zi_grid, d.se_b, d.ke_b, d.wv_b, d.wl_b, d.se_a, d.ke_a, d.wv_a, d.wl_a, d.wthl_sfc, d.wqw_sfc, d.rho_zt, d.tke, d.pint, d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_integrals(ShocEnergyIntegralsData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_integrals_c(d.shcol, d.nlev, d.host_dse, d.pdel, d.rtm, d.rcm, d.u_wind, d.v_wind, d.se_int, d.ke_int, d.wv_int, d.wl_int);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_total_fixer(ShocEnergyTotalFixerData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_total_fixer_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.zt_grid, d.zi_grid, d.se_b, d.ke_b, d.wv_b, d.wl_b, d.se_a, d.ke_a, d.wv_a, d.wl_a, d.wthl_sfc, d.wqw_sfc, d.rho_zt, d.te_a, d.te_b);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_threshold_fixer(ShocEnergyThresholdFixerData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_threshold_fixer_c(d.shcol, d.nlev, d.nlevi, d.pint, d.tke, d.te_a, d.te_b, d.se_dis, d.shoctop);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_dse_fixer(ShocEnergyDseFixerData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_dse_fixer_c(d.shcol, d.nlev, d.se_dis, d.shoctop, d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void calc_shoc_vertflux(CalcShocVertfluxData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar, d.vertflux);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void calc_shoc_varorcovar(CalcShocVarorcovarData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi, d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_tmpi(ComputeTmpiData& d)
{
  shoc_init(d.nlevi - 1, true); // nlev = nlevi - 1
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_tmpi_c(d.nlevi, d.shcol, d.dtime, d.rho_zi, d.dz_zi, d.tmpi);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void dp_inverse(DpInverseData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  dp_inverse_c(d.nlev, d.shcol, d.rho_zt, d.dz_zt, d.rdp_zt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void sfc_fluxes(SfcFluxesData& d)
{
  shoc_init(1, true); // single layer function
  d.transpose<ekat::TransposeDirection::c2f>();
  sfc_fluxes_c(d.shcol, d.num_tracer, d.dtime, d.rho_zi_sfc, d.rdp_zt_sfc, d.wthl_sfc, d.wqw_sfc, d.wtke_sfc, d.wtracer_sfc, d.thetal, d.qw, d.tke, d.wtracer);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void impli_srf_stress_term(ImpliSrfStressTermData& d)
{
  shoc_init(1, true); // single layer function
  impli_srf_stress_term_c(d.shcol, d.rho_zi_sfc, d.uw_sfc, d.vw_sfc, d.u_wind_sfc, d.v_wind_sfc, d.ksrf);
}

void tke_srf_flux_term(TkeSrfFluxTermData& d)
{
  shoc_init(1, true); // single layer function
  tke_srf_flux_term_c(d.shcol, d.uw_sfc, d.vw_sfc, d.wtke_sfc);
}

void integ_column_stability(IntegColumnStabilityData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  integ_column_stability_c(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void check_tke(CheckTkeData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  check_tke_c(d.shcol, d.nlev, d.tke);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_tke(ShocTkeData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_tke_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.wthv_sec, d.shoc_mix, d.dz_zi, d.dz_zt, d.pres, d.u_wind, d.v_wind, d.brunt, d.obklen, d.zt_grid, d.zi_grid, d.pblh, d.tke, d.tk, d.tkh, d.isotropy);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_shr_prod(ComputeShrProdData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_shr_prod_c(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind, d.v_wind, d.sterm);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void isotropic_ts(IsotropicTsData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  isotropic_ts_c(d.nlev, d.shcol, d.brunt_int, d.tke, d.a_diss, d.brunt, d.isotropy);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void adv_sgs_tke(AdvSgsTkeData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  adv_sgs_tke_c(d.nlev, d.shcol, d.dtime, d.shoc_mix, d.wthv_sec, d.sterm_zt, d.tk, d.tke, d.a_diss);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void eddy_diffusivities(EddyDiffusivitiesData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  eddy_diffusivities_c(d.nlev, d.shcol, d.obklen, d.pblh, d.zt_grid, d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_length(ShocLengthData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_length_c(d.shcol, d.nlev, d.nlevi, d.host_dx, d.host_dy, d.zt_grid, d.zi_grid, d.dz_zt, d.tke, d.thv, d.brunt, d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_brunt_shoc_length(ComputeBruntShocLengthData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_brunt_shoc_length_c(d.nlev, d.nlevi, d.shcol, d.dz_zt, d.thv, d.thv_zi, d.brunt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_l_inf_shoc_length(ComputeLInfShocLengthData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_l_inf_shoc_length_c(d.nlev, d.shcol, d.zt_grid, d.dz_zt, d.tke, d.l_inf);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_shoc_mix_shoc_length(ComputeShocMixShocLengthData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_shoc_mix_shoc_length_c(d.nlev, d.shcol, d.tke, d.brunt, d.zt_grid, d.l_inf, d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void check_length_scale_shoc_length(CheckLengthScaleShocLengthData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  check_length_scale_shoc_length_c(d.nlev, d.shcol, d.host_dx, d.host_dy, d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void fterms_input_for_diag_third_shoc_moment(FtermsInputForDiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  fterms_input_for_diag_third_shoc_moment_c(d.dz_zi, d.dz_zt, d.dz_zt_kc, d.isotropy_zi, d.brunt_zi, d.thetal_zi, &d.thedz, &d.thedz2, &d.iso, &d.isosqrd, &d.buoy_sgs2, &d.bet2);
}

void aa_terms_diag_third_shoc_moment(AaTermsDiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  aa_terms_diag_third_shoc_moment_c(d.omega0, d.omega1, d.omega2, d.x0, d.x1, d.y0, d.y1, &d.aa0, &d.aa1);
}

void f0_to_f5_diag_third_shoc_moment(F0ToF5DiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  f0_to_f5_diag_third_shoc_moment_c(d.thedz, d.thedz2, d.bet2, d.iso, d.isosqrd, d.wthl_sec, d.wthl_sec_kc, d.wthl_sec_kb, d.thl_sec_kc, d.thl_sec_kb, d.w_sec, d.w_sec_kc, d.w_sec_zi, d.tke, d.tke_kc, &d.f0, &d.f1, &d.f2, &d.f3, &d.f4, &d.f5);
}

void omega_terms_diag_third_shoc_moment(OmegaTermsDiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  omega_terms_diag_third_shoc_moment_c(d.buoy_sgs2, d.f3, d.f4, &d.omega0, &d.omega1, &d.omega2);
}

void x_y_terms_diag_third_shoc_moment(XYTermsDiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  x_y_terms_diag_third_shoc_moment_c(d.buoy_sgs2, d.f0, d.f1, d.f2, &d.x0, &d.y0, &d.x1, &d.y1);
}

void w3_diag_third_shoc_moment(W3DiagThirdShocMomentData& d)
{
  shoc_init(1, true); // single level function
  w3_diag_third_shoc_moment_c(d.aa0, d.aa1, d.x0, d.x1, d.f5, &d.w3);
}

void clipping_diag_third_shoc_moments(ClippingDiagThirdShocMomentsData& d)
{
  shoc_init(d.nlevi - 1, true); // nlev = nlevi - 1
  d.transpose<ekat::TransposeDirection::c2f>();
  clipping_diag_third_shoc_moments_c(d.nlevi, d.shcol, d.w_sec_zi, d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_second_moments_srf(DiagSecondMomentsSrfData& d)
{
  shoc_init(1, true); // single level function
  shoc_diag_second_moments_srf_c(d.shcol, d.wthl_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar);
}

void linear_interp(LinearInterpData& d)
{
  shoc_init(d.km1, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  linear_interp_c(d.x1, d.x2, d.y1, d.y2, d.km1, d.km2, d.ncol, d.minthresh);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_third_shoc_moments(DiagThirdShocMomentsData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_third_shoc_moments_c(d.shcol, d.nlev, d.nlevi, d.w_sec, d.thl_sec, d.wthl_sec, d.isotropy, d.brunt, d.thetal, d.tke, d.dz_zt, d.dz_zi, d.zt_grid, d.zi_grid, d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_diag_third_shoc_moment(ComputeDiagThirdShocMomentData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_diag_third_shoc_moment_c(d.shcol, d.nlev, d.nlevi, d.w_sec, d.thl_sec, d.wthl_sec, d.tke, d.dz_zt, d.dz_zi, d.isotropy_zi, d.brunt_zi, d.w_sec_zi, d.thetal_zi, d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_assumed_pdf(ShocAssumedPdfData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_assumed_pdf_c(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.w_field, d.thl_sec, d.qw_sec, d.wthl_sec, d.w_sec, d.wqw_sec, d.qwthl_sec, d.w3, d.pres, d.zt_grid, d.zi_grid, d.shoc_cldfrac, d.shoc_ql, d.wqls, d.wthv_sec, d.shoc_ql2);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_assumed_pdf_tilde_to_real(ShocAssumedPdfTildeToRealData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_tilde_to_real_c(d.w_first, d.sqrtw2, &d.w1);
}

void shoc_assumed_pdf_vv_parameters(ShocAssumedPdfVvParametersData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_vv_parameters_c(d.w_first, d.w_sec, d.w3var, &d.skew_w, &d.w1_1, &d.w1_2, &d.w2_1, &d.w2_2, &d.a);
}

void shoc_assumed_pdf_thl_parameters(ShocAssumedPdfThlParametersData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_thl_parameters_c(d.wthlsec, d.sqrtw2, d.sqrtthl, d.thlsec, d.thl_first, d.w1_1, d.w1_2, d.skew_w, d.a, d.dothetal_skew, &d.thl1_1, &d.thl1_2, &d.thl2_1, &d.thl2_2, &d.sqrtthl2_1, &d.sqrtthl2_2);
}

void shoc_assumed_pdf_qw_parameters(ShocAssumedPdfQwParametersData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_qw_parameters_c(d.wqwsec, d.sqrtw2, d.skew_w, d.sqrtqt, d.qwsec, d.w1_2, d.w1_1, d.qw_first, d.a, &d.qw1_1, &d.qw1_2, &d.qw2_1, &d.qw2_2, &d.sqrtqw2_1, &d.sqrtqw2_2);
}

void shoc_assumed_pdf_inplume_correlations(ShocAssumedPdfInplumeCorrelationsData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_inplume_correlations_c(d.sqrtqw2_1, d.sqrtthl2_1, d.a, d.sqrtqw2_2, d.sqrtthl2_2, d.qwthlsec, d.qw1_1, d.qw_first, d.thl1_1, d.thl_first, d.qw1_2, d.thl1_2, &d.r_qwthl_1);
}

void shoc_assumed_pdf_compute_temperature(ShocAssumedPdfComputeTemperatureData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_temperature_c(d.thl1, d.basepres, d.pval, &d.tl1);
}

void shoc_assumed_pdf_compute_qs(ShocAssumedPdfComputeQsData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_qs_c(d.tl1_1, d.tl1_2, d.pval, &d.qs1, &d.beta1, &d.qs2, &d.beta2);
}

void shoc_assumed_pdf_compute_s(ShocAssumedPdfComputeSData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_s_c(d.qw1, d.qs1, d.beta, d.pval, d.thl2, d.qw2, d.sqrtthl2, d.sqrtqw2, d.r_qwthl, &d.s, &d.std_s, &d.qn, &d.c);
}

void shoc_assumed_pdf_compute_sgs_liquid(ShocAssumedPdfComputeSgsLiquidData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_sgs_liquid_c(d.a, d.ql1, d.ql2, &d.shoc_ql);
}

void shoc_assumed_pdf_compute_cloud_liquid_variance(ShocAssumedPdfComputeCloudLiquidVarianceData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_cloud_liquid_variance_c(d.a, d.s1, d.ql1, d.c1, d.std_s1, d.s2, d.ql2, d.c2, d.std_s2, d.shoc_ql, &d.shoc_ql2);
}

void shoc_assumed_pdf_compute_liquid_water_flux(ShocAssumedPdfComputeLiquidWaterFluxData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_liquid_water_flux_c(d.a, d.w1_1, d.w_first, d.ql1, d.w1_2, d.ql2, &d.wqls);
}

void shoc_assumed_pdf_compute_buoyancy_flux(ShocAssumedPdfComputeBuoyancyFluxData& d)
{
  shoc_init(1, true); // single level function
  shoc_assumed_pdf_compute_buoyancy_flux_c(d.wthlsec, d.epsterm, d.wqwsec, d.pval, d.wqls, &d.wthv_sec);
}

void diag_second_moments_ubycond(DiagSecondMomentsUbycondData& d)
{
  shoc_init(1, true); // single level function
  shoc_diag_second_moments_ubycond_c(d.shcol, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec);
}

void pblintd_init_pot(PblintdInitPotData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_pblintd_init_pot_c(d.shcol, d.nlev, d.thl, d.ql, d.q, d.thv);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void pblintd_cldcheck(PblintdCldcheckData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_pblintd_cldcheck_c(d.shcol, d.nlev, d.nlevi, d.zi, d.cldn, d.pblh);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_second_moments_lbycond(DiagSecondMomentsLbycondData& d)
{
  shoc_init(1, true); // single level function
  diag_second_moments_lbycond_c(d.shcol, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar, d.wthl_sec, d.wqw_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.thl_sec, d.qw_sec, d.qwthl_sec);
}

void diag_second_moments(DiagSecondMomentsData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_second_moments_c(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh, d.tk,
       d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec,
       d.vw_sec, d.wtke_sec, d.w_sec);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_second_shoc_moments(DiagSecondShocMomentsData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_second_shoc_moments_c(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh,
      d.tk, d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.thl_sec, d.qw_sec,
      d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.w_sec);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_shoc_vapor(ComputeShocVaporData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_shoc_vapor_c(d.shcol, d.nlev, d.qw, d.ql, d.qv);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void update_prognostics_implicit(UpdatePrognosticsImplicitData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  update_prognostics_implicit_c(d.shcol, d.nlev, d.nlevi, d.num_tracer, d.dtime,
                                d.dz_zt, d.dz_zi, d.rho_zt, d.zt_grid, d.zi_grid,
                                d.tk, d.tkh, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc,
                                d.wtracer_sfc, d.thetal, d.qw, d.tracer, d.tke, d.u_wind, d.v_wind);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_main(ShocMainData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_main_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.host_dx, d.host_dy, d.thv, d.zt_grid, d.zi_grid,
              d.pres, d.presi, d.pdel, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.wtracer_sfc,
              d.num_qtracers, d.w_field, d.inv_exner, d.phis, d.host_dse, d.tke, d.thetal, d.qw,
              d.u_wind, d.v_wind, d.qtracers, d.wthv_sec, d.tkh, d.tk, d.shoc_ql, d.shoc_cldfrac, d.pblh,
              d.shoc_mix, d.isotropy, d.w_sec, d.thl_sec, d.qw_sec, d.qwthl_sec, d.wthl_sec, d.wqw_sec,
              d.wtke_sec, d.uw_sec, d.vw_sec, d.w3, d.wqls_sec, d.brunt, d.shoc_ql2, &d.elapsed_s);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_main_with_init(ShocMainData& d)
{
  using C = scream::physics::Constants<Real>;

  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_init_for_main_bfb_c(d.nlev, C::gravit, C::Rair, C::RH2O, C::Cpair, C::ZVIR, C::LatVap, C::LatIce, C::Karman,
                           d.pref_mid, d.nbot_shoc, d.ntop_shoc+1);
  shoc_use_cxx_c(false);


  shoc_main_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.host_dx, d.host_dy, d.thv, d.zt_grid, d.zi_grid,
              d.pres, d.presi, d.pdel, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.wtracer_sfc, d.num_qtracers,
              d.w_field, d.inv_exner, d.phis, d.host_dse, d.tke, d.thetal, d.qw, d.u_wind, d.v_wind, d.qtracers,
              d.wthv_sec, d.tkh, d.tk, d.shoc_ql, d.shoc_cldfrac, d.pblh, d.shoc_mix, d.isotropy, d.w_sec,
              d.thl_sec, d.qw_sec, d.qwthl_sec, d.wthl_sec, d.wqw_sec, d.wtke_sec, d.uw_sec, d.vw_sec, d.w3,
              d.wqls_sec, d.brunt, d.shoc_ql2, &d.elapsed_s);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void pblintd_height(PblintdHeightData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  pblintd_height_c(d.shcol, d.nlev, d.npbl, d.z, d.u, d.v, d.ustar, d.thv, d.thv_ref, d.pblh, d.rino, d.check);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void vd_shoc_decomp_and_solve(VdShocDecompandSolveData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  // Call decomp subroutine
  vd_shoc_decomp_c(d.shcol, d.nlev, d.nlevi, d.kv_term, d.tmpi, d.rdp_zt, d.dtime, d.flux, d.du, d.dl, d.d);
  // Call solver for each problem. The `var` array represents 3d
  // data with an entry per (shcol, nlev, n_rhs). Fortran requires
  // 2d data (shcol, nlev) for each rhs.
  const Int size = d.shcol*d.nlev;
  for (Int n=0; n<d.n_rhs; ++n) {
    // Copy var to rhs
    for(Int s=0; s<size; ++s) {
      d.rhs[s] = d.var[n*size+s];
    }
    vd_shoc_solve_c(d.shcol, d.nlev, d.du, d.dl, d.d, d.rhs);
    // Copy rhs to var
    for(Int s=0; s<size; ++s) {
      d.var[n*size+s] = d.rhs[s];
    }
  }
  d.transpose<ekat::TransposeDirection::f2c>();
}

void pblintd_surf_temp(PblintdSurfTempData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  pblintd_surf_temp_c(d.shcol, d.nlev, d.nlevi, d.z, d.ustar, d.obklen, d.kbfs, d.thv, d.tlv, d.pblh, d.check, d.rino);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void pblintd_check_pblh(PblintdCheckPblhData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  pblintd_check_pblh_c(d.shcol, d.nlev, d.nlevi, d.z, d.ustar, d.check, d.pblh);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void pblintd(PblintdData& d)
{
  shoc_init(d.nlev, true, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  pblintd_c(d.shcol, d.nlev, d.nlevi, d.npbl, d.z, d.zi, d.thl, d.ql, d.q, d.u, d.v, d.ustar, d.obklen, d.kbfs, d.cldn, d.pblh);
  d.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

//
// _f function definitions. These expect data in C layout
//

void calc_shoc_varorcovar_f(Int shcol, Int nlev, Int nlevi, Real tunefac,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  ekat::device_to_host({varorcovar}, shcol, nlevi, inout_views, true);
}

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  ekat::device_to_host({vertflux}, shcol, nlevi, inout_views, true);
}

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar)
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

void shoc_diag_second_moments_ubycond_f(Int shcol, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec,
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

void update_host_dse_f(Int shcol, Int nlev, Real* thlm, Real* shoc_ql, Real* inv_exner, Real* zt_grid,
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
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d, true);

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
  ekat::device_to_host({host_dse}, shcol, nlev, inout_views, true);
}

void compute_diag_third_shoc_moment_f(Int shcol, Int nlev, Int nlevi, Real* w_sec,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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

    SHF::compute_diag_third_shoc_moment(team, nlev, nlevi, w_sec_s, thl_sec_s,
                                        wthl_sec_s, tke_s, dz_zt_s, dz_zi_s, isotropy_zi_s,
                                        brunt_zi_s, w_sec_zi_s, thetal_zi_s, w3_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {w3_d};
  ekat::device_to_host({w3}, shcol, nlevi, inout_views, true);
}

void shoc_pblintd_init_pot_f(Int shcol, Int nlev, Real *thl, Real* ql, Real* q,
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
  ekat::host_to_device({thl, ql, q}, shcol, nlev, temp_d, true);

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
  ekat::device_to_host({thv}, shcol, nlev, inout_views, true);
}

void compute_shoc_mix_shoc_length_f(Int nlev, Int shcol, Real* tke, Real* brunt,
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
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d, true);

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

    SHF::compute_shoc_mix_shoc_length(team, nlev, tke_s, brunt_s, zt_grid_s, l_inf_s,
                                      shoc_mix_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {shoc_mix_d};
  ekat::device_to_host({shoc_mix}, shcol, nlev, inout_views, true);
}

void check_tke_f(Int shcol, Int nlev, Real* tke)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  std::vector<view_2d> temp_2d_d(1);

  // Sync to device
  ekat::host_to_device({tke}, shcol, nlev, temp_2d_d, true);

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
  ekat::device_to_host({tke}, shcol, nlev, inout_views, true);
}

void linear_interp_f(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh)
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

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
  ekat::device_to_host({y2}, ncol, km2, inout_views, true);
}

void clipping_diag_third_shoc_moments_f(Int nlevi, Int shcol, Real *w_sec_zi,
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
  ekat::host_to_device({w_sec_zi, w3}, shcol, nlevi, temp_d, true);

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
  ekat::device_to_host({w3}, shcol, nlevi, inout_views, true);
}

void shoc_energy_integrals_f(Int shcol, Int nlev, Real *host_dse, Real *pdel,
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
  ekat::host_to_device(ptr_array, shcol, nlev, temp_d, true);

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

void diag_second_moments_lbycond_f(Int shcol, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar,
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

void diag_second_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind,
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

  ekat::host_to_device(ptr_array, dim1_array, dim2_array, temp_2d, true);

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

    SHOC::diag_second_moments(team, nlev, nlevi, thetal_1d, qw_1d, u_wind_1d, v_wind_1d, tke_1d, isotropy_1d, tkh_1d, tk_1d,
                     dz_zi_1d, zt_grid_1d, zi_grid_1d, shoc_mix_1d, isotropy_zi_1d, tkh_zi_1d, tk_zi_1d,
                     thl_sec_1d, qw_sec_1d, wthl_sec_1d, wqw_sec_1d,
                     qwthl_sec_1d, uw_sec_1d, vw_sec_1d, wtke_sec_1d, w_sec_1d);


  });

  std::vector<Int> dim1(9, shcol);
  std::vector<Int> dim2 = {nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlev };
  std::vector<view_2d> host_views = {thl_sec_2d, qw_sec_2d, wthl_sec_2d, wqw_sec_2d, qwthl_sec_2d, uw_sec_2d, vw_sec_2d, wtke_sec_2d, w_sec_2d};
  ekat::device_to_host({thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec}, dim1, dim2, host_views, true);
}

void diag_second_shoc_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke,
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

  ekat::host_to_device(ptr_array, dim1_array, dim2_array, temp_2d, true);

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
       thetal_1d, qw_1d, u_wind_1d, v_wind_1d, tke_1d, isotropy_1d, tkh_1d, tk_1d, dz_zi_1d, zt_grid_1d, zi_grid_1d, shoc_mix_1d,
       wthl_s, wqw_s, uw_s, vw_s, ustar2_s, wstar_s,
       workspace, thl_sec_1d, qw_sec_1d, wthl_sec_1d, wqw_sec_1d, qwthl_sec_1d,
       uw_sec_1d, vw_sec_1d, wtke_sec_1d, w_sec_1d);
  });

  std::vector<Int> dim1(9, shcol);
  std::vector<Int> dim2 = {nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlevi, nlev };
  std::vector<view_2d> host_2d_views = {thl_sec_2d, qw_sec_2d, wthl_sec_2d, wqw_sec_2d, qwthl_sec_2d, uw_sec_2d, vw_sec_2d, wtke_sec_2d, w_sec_2d};
  ekat::device_to_host({thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec}, dim1, dim2, host_2d_views, true);
}

void compute_brunt_shoc_length_f(Int nlev, Int nlevi, Int shcol, Real* dz_zt, Real* thv, Real* thv_zi, Real* brunt)
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  ekat::device_to_host({brunt}, shcol, nlev, inout_views, true);
}

void compute_l_inf_shoc_length_f(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
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
  ekat::host_to_device({zt_grid, dz_zt, tke}, shcol, nlev, temp_d, true);

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

void check_length_scale_shoc_length_f(Int nlev, Int shcol, Real* host_dx, Real* host_dy, Real* shoc_mix)
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
  ekat::host_to_device({shoc_mix}, shcol, nlev, temp_2d_d, true);

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
  ekat::device_to_host({shoc_mix}, shcol, nlev, inout_views, true);
}

void shoc_diag_obklen_f(Int shcol, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc, Real* thl_sfc,
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

void shoc_pblintd_cldcheck_f(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh) {
  using SHOC    = Functions<Real, DefaultDevice>;
  using Spack   = typename SHOC::Spack;
  using Scalar  = typename SHOC::Scalar;
  using view_2d = typename SHOC::view_2d<Spack>;
  using view_1d = typename SHOC::view_1d<Scalar>;

  std::vector<Int> dim1(2, shcol);
  std::vector<Int> dim2 = {nlevi, nlev};

  std::vector<view_2d> cldcheck_2d(2);
  ekat::host_to_device({zi, cldn}, dim1, dim2, cldcheck_2d, true);

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

void shoc_length_f(Int shcol, Int nlev, Int nlevi, Real* host_dx, Real* host_dy,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

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

    SHF::shoc_length(team,nlev,nlevi,host_dx_s,host_dy_s,
                     zt_grid_s,zi_grid_s,dz_zt_s,tke_s,
                     thv_s,workspace,brunt_s,shoc_mix_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {brunt_d,shoc_mix_d};
  ekat::device_to_host({brunt,shoc_mix}, shcol, nlev, inout_views, true);
}

void shoc_energy_fixer_f(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* zt_grid,
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
  ekat::host_to_device(ptr_array_2d, dim1_sizes, dim2_sizes, temp_2d_d, true);

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
  ekat::device_to_host({host_dse}, shcol, nlev, inout_views, true);
}

void compute_shoc_vapor_f(Int shcol, Int nlev, Real* qw, Real* ql, Real* qv)
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
  ekat::host_to_device( {qw,  ql, qv}, shcol, nlev, temp_d, true);

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
  ekat::device_to_host({qv}, shcol, nlev, inout_views, true);
}

void update_prognostics_implicit_f(Int shcol, Int nlev, Int nlevi, Int num_tracer, Real dtime,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);
  ekat::host_to_device({tracer}, shcol, nlev, num_tracer, temp_3d_d, true);

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
  ekat::device_to_host({thetal, qw, u_wind, v_wind, tke}, shcol, nlev, inout_views_2d, true);

  std::vector<view_3d> inout_views = {qtracers_f90_d};
  ekat::device_to_host({tracer}, shcol, nlev, num_tracer, inout_views, true);
}

void diag_third_shoc_moments_f(Int shcol, Int nlev, Int nlevi, Real* w_sec, Real* thl_sec,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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

    SHF::diag_third_shoc_moments(team, nlev, nlevi, wsec_s, thl_sec_s,
                                 wthl_sec_s, isotropy_s, brunt_s, thetal_s, tke_s,
                                 dz_zt_s, dz_zi_s, zt_grid_s, zi_grid_s,
                                 workspace,
                                 w3_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {w3_d};
  ekat::device_to_host({w3}, shcol, nlevi, inout_views, true);
}

void adv_sgs_tke_f(Int nlev, Int shcol, Real dtime, Real* shoc_mix, Real* wthv_sec,
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
  ekat::host_to_device(ptr_array, shcol, nlev, temp_d, true);

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
  ekat::device_to_host({tke, a_diss}, shcol, nlev, inout_views, true);
}

void shoc_assumed_pdf_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* w_field,
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  ekat::device_to_host({shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2}, shcol, nlev, out_views, true);
}
void compute_shr_prod_f(Int nlevi, Int nlev, Int shcol, Real* dz_zi, Real* u_wind, Real* v_wind, Real* sterm)
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  ekat::device_to_host({sterm}, shcol, nlevi, inout_views, true);
}

void compute_tmpi_f(Int nlevi, Int shcol, Real dtime, Real *rho_zi, Real *dz_zi, Real *tmpi)
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
  ekat::host_to_device({rho_zi,  dz_zi, tmpi}, shcol, nlevi, temp_d, true);

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
  ekat::device_to_host({tmpi}, shcol, nlevi, inout_views, true);
}

void integ_column_stability_f(Int nlev, Int shcol, Real *dz_zt,
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
  ekat::host_to_device({dz_zt, pres, brunt}, shcol, nlev, temp_d, true);

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

void isotropic_ts_f(Int nlev, Int shcol, Real* brunt_int, Real* tke,
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
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d, true);

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

      SHF::isotropic_ts(team, nlev, brunt_int_s, tke_s, a_diss_s, brunt_s, isotropy_s);
    });

  // Sync back to host
  std::vector<view_2d> inout_views = {isotropy_d};
  ekat::device_to_host({isotropy}, shcol, nlev, inout_views, true);

}

void dp_inverse_f(Int nlev, Int shcol, Real *rho_zt, Real *dz_zt, Real *rdp_zt)
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
  ekat::host_to_device({rho_zt,  dz_zt, rdp_zt}, shcol, nlev, temp_d, true);

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
  ekat::device_to_host({rdp_zt}, shcol, nlev, inout_views, true);
}

int shoc_init_f(Int nlev, Real *pref_mid, Int nbot_shoc, Int ntop_shoc)
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

Int shoc_main_f(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Int npbl, Real* host_dx, Real* host_dy, Real* thv, Real* zt_grid,
                Real* zi_grid, Real* pres, Real* presi, Real* pdel, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc,
                Real* wtracer_sfc, Int num_qtracers, Real* w_field, Real* inv_exner, Real* phis, Real* host_dse, Real* tke,
                Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk,
                Real* shoc_ql, Real* shoc_cldfrac, Real* pblh, Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec,
                Real* qw_sec, Real* qwthl_sec, Real* wthl_sec, Real* wqw_sec, Real* wtke_sec, Real* uw_sec, Real* vw_sec,
                Real* w3, Real* wqls_sec, Real* brunt, Real* shoc_ql2)
{
  // tkh is a local variable in C++ impl
  (void)tkh;

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
  static constexpr Int num_2d_arrays = 34;
  static constexpr Int num_3d_arrays = 1;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);
  std::vector<view_3d> temp_3d_d(num_3d_arrays);

  std::vector<int> dim1_2d_sizes = {shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol,
                                    shcol, shcol, shcol, shcol, shcol};
  std::vector<int> dim2_2d_sizes = {nlev,  nlevi, nlev,         nlevi, nlev,
                                    nlev,  nlev,  num_qtracers, nlev,  nlev,
                                    nlev,  nlev,  nlev,         nlev,  nlev,
                                    nlev,  nlev,  nlev,  nlev,
                                    nlev,  nlev,  nlev,         nlevi, nlevi,
                                    nlevi, nlevi, nlevi,        nlevi, nlevi,
                                    nlevi, nlevi, nlev,         nlev,  nlev};

  std::vector<const Real*> ptr_array_1d = {host_dx, host_dy, wthl_sfc, wqw_sfc,
                                           uw_sfc,  vw_sfc,  phis};
  std::vector<const Real*> ptr_array_2d = {zt_grid,   zi_grid,  pres,        presi,        pdel,
                                           thv,       w_field,  wtracer_sfc, inv_exner,        host_dse,
                                           tke,       thetal,   qw,          u_wind,       v_wind,
                                           wthv_sec,  tk,       shoc_cldfrac, shoc_ql,
                                           shoc_ql2,  shoc_mix, w_sec,       thl_sec,      qw_sec,
                                           qwthl_sec, wthl_sec, wqw_sec,     wtke_sec,     uw_sec,
                                           vw_sec,    w3,       wqls_sec,    brunt,        isotropy};

  ScreamDeepCopy::copy_to_device(ptr_array_1d, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array_2d, dim1_2d_sizes, dim2_2d_sizes, temp_2d_d, true);
  ekat::host_to_device({qtracers}, shcol, nlev, num_qtracers, temp_3d_d, true);

  Int index_counter = 0;
  view_1d
    host_dx_d (temp_1d_d[index_counter++]),
    host_dy_d (temp_1d_d[index_counter++]),
    wthl_sfc_d(temp_1d_d[index_counter++]),
    wqw_sfc_d (temp_1d_d[index_counter++]),
    uw_sfc_d  (temp_1d_d[index_counter++]),
    vw_sfc_d  (temp_1d_d[index_counter++]),
    phis_d    (temp_1d_d[index_counter++]),
    pblh_d    ("pblh",shcol);

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
  SHF::SHOCOutput shoc_output{pblh_d, shoc_ql2_d};
  SHF::SHOCHistoryOutput shoc_history_output{shoc_mix_d,  w_sec_d,    thl_sec_d, qw_sec_d,
                                             qwthl_sec_d, wthl_sec_d, wqw_sec_d, wtke_sec_d,
                                             uw_sec_d,    vw_sec_d,   w3_d,      wqls_sec_d,
                                             brunt_d,     isotropy_d};

  const auto nlevi_packs = ekat::npack<Spack>(nlevi);

#ifdef SCREAM_SMALL_KERNELS
  view_1d
    se_b   ("se_b", shcol),
    ke_b   ("ke_b", shcol),
    wv_b   ("wv_b", shcol),
    wl_b   ("wl_b", shcol),
    se_a   ("se_a", shcol),
    ke_a   ("ke_a", shcol),
    wv_a   ("wv_a", shcol),
    wl_a   ("wl_a", shcol),
    ustar  ("ustar", shcol),
    kbfs   ("kbfs", shcol),
    obklen ("obklen", shcol),
    ustar2 ("ustar2", shcol),
    wstar  ("wstar", shcol);

  view_2d
    rho_zt  ("rho_zt",  shcol, nlevi_packs),
    shoc_qv ("shoc_qv", shcol, nlevi_packs),
    dz_zt   ("dz_zt",   shcol, nlevi_packs),
    dz_zi   ("dz_zi",   shcol, nlevi_packs),
    tkhv    ("tkh",     shcol, nlevi_packs);

  SHF::SHOCTemporaries shoc_temporaries{
    se_b, ke_b, wv_b, wl_b, se_a, ke_a, wv_a, wl_a, ustar, kbfs, obklen, ustar2, wstar,
    rho_zt, shoc_qv, dz_zt, dz_zi, tkhv};
#endif

  // Create local workspace
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
  ekat::WorkspaceManager<Spack, SHF::KT::Device> workspace_mgr(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);

  const auto elapsed_microsec = SHF::shoc_main(shcol, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
                                               workspace_mgr,
                                               shoc_input, shoc_input_output, shoc_output, shoc_history_output
#ifdef SCREAM_SMALL_KERNELS
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
  ekat::device_to_host(ptr_array_2d_out, dim1_2d_out, dim2_2d_out, out_views_2d, true);

  // 3d
  std::vector<view_3d> out_views_3d = {qtracers_f90_d};
  ekat::device_to_host({qtracers}, shcol, nlev, num_qtracers, out_views_3d, true);

  return elapsed_microsec;
}

void pblintd_height_f(Int shcol, Int nlev, Int npbl, Real* z, Real* u, Real* v, Real* ustar, Real* thv, Real* thv_ref, Real* pblh, Real* rino, bool* check)
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
  ekat::host_to_device({z, u, v, thv, rino}, shcol, nlev, views_2d, true);

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
  ekat::device_to_host({rino}, shcol, nlev, out_2d_views, true);

  std::vector<bview_1d> out_bool_1d_views = {check_1d};
  ScreamDeepCopy::copy_to_host({check}, shcol, out_bool_1d_views);
}

void vd_shoc_decomp_and_solve_f(Int shcol, Int nlev, Int nlevi, Int num_rhs, Real* kv_term, Real* tmpi, Real* rdp_zt, Real dtime,
                                Real* flux, Real* var)
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);
  ekat::host_to_device({var}, shcol, nlev, num_rhs, temp_3d_d, true);

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
  ekat::device_to_host({var}, shcol, nlev, num_rhs, inout_views, true);
}

void shoc_grid_f(Int shcol, Int nlev, Int nlevi, Real* zt_grid, Real* zi_grid, Real* pdel, Real* dz_zt, Real* dz_zi, Real* rho_zt)
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
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

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
  ekat::device_to_host<Int>({dz_zt, dz_zi, rho_zt}, {shcol, shcol, shcol}, {nlev, nlevi, nlev}, inout_views, true);
}

void eddy_diffusivities_f(Int nlev, Int shcol, Real* obklen, Real* pblh, Real* zt_grid, Real* shoc_mix, Real* sterm_zt,
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

  static constexpr Int num_1d_arrays = 2;
  static constexpr Int num_2d_arrays = 7;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);

  std::vector<const Real*> ptr_array = {zt_grid, shoc_mix, sterm_zt, isotropy,
                                        tke,     tkh,      tk};

  // Sync to device
  ScreamDeepCopy::copy_to_device({obklen, pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, shcol, nlev, temp_2d_d, true);

  view_1d
    obklen_d(temp_1d_d[0]),
    pblh_d(temp_1d_d[1]);

  view_2d
    zt_grid_d(temp_2d_d[0]),
    shoc_mix_d(temp_2d_d[1]),
    sterm_zt_d(temp_2d_d[2]),
    isotropy_d(temp_2d_d[3]),
    tke_d(temp_2d_d[4]),
    tkh_d(temp_2d_d[5]),
    tk_d(temp_2d_d[6]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar obklen_s{obklen_d(i)};
    const Scalar pblh_s{pblh_d(i)};

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);
    const auto sterm_zt_s = ekat::subview(sterm_zt_d, i);
    const auto isotropy_s = ekat::subview(isotropy_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto tkh_s = ekat::subview(tkh_d, i);
    const auto tk_s = ekat::subview(tk_d, i);

    SHF::eddy_diffusivities(team, nlev, obklen_s, pblh_s, zt_grid_s, shoc_mix_s, sterm_zt_s, isotropy_s, tke_s, tkh_s, tk_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tkh_d, tk_d};
  ekat::device_to_host({tkh, tk}, shcol, nlev, inout_views, true);
}

void pblintd_surf_temp_f(Int shcol, Int nlev, Int nlevi, Real* z, Real* ustar, Real* obklen, Real* kbfs, Real* thv, Real* tlv, Real* pblh, bool* check, Real* rino)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using Scalar     = typename SHOC::Scalar;
  using view_1d      = typename SHOC::view_1d<Scalar>;
  using view_bool_1d = typename SHOC::view_1d<bool>;
  using view_2d      = typename SHOC::view_2d<Spack>;

  std::vector<view_2d> views_2d(3);
  ekat::host_to_device({z, thv, rino}, shcol, nlev, views_2d, true);
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
  ekat::device_to_host({rino}, shcol, nlev, out_2d_views, true);

  std::vector<view_bool_1d> out_bool_1d_views = {check_1d};
  ScreamDeepCopy::copy_to_host({check}, shcol, out_bool_1d_views);
}

void pblintd_check_pblh_f(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* ustar, bool* check, Real* pblh)
{
  using SHOC         = Functions<Real, DefaultDevice>;
  using Spack        = typename SHOC::Spack;
  using Scalar       = typename SHOC::Scalar;
  using view_bool_1d = typename SHOC::view_1d<bool>;
  using view_1d      = typename SHOC::view_1d<Scalar>;
  using view_2d      = typename SHOC::view_2d<Spack>;

  std::vector<view_2d> views_2d(1);
  ekat::host_to_device({z}, shcol, nlev, views_2d, true);
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

void pblintd_f(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* zi, Real* thl, Real* ql, Real* q, Real* u, Real* v, Real* ustar, Real* obklen, Real* kbfs, Real* cldn, Real* pblh)
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
    ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

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

void shoc_tke_f(Int shcol, Int nlev, Int nlevi, Real dtime, Real* wthv_sec, Real* shoc_mix, Real* dz_zi, Real* dz_zt, Real* pres,
                Real* u_wind, Real* v_wind, Real* brunt, Real* obklen, Real* zt_grid, Real* zi_grid, Real* pblh, Real* tke, Real* tk,
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

  static constexpr Int num_1d_arrays = 2;
  static constexpr Int num_2d_arrays = 14;

  std::vector<view_1d> temp_1d_d(num_1d_arrays);
  std::vector<view_2d> temp_2d_d(num_2d_arrays);

  std::vector<int> dim1_sizes(num_2d_arrays, shcol);
  std::vector<int> dim2_sizes = {nlev, nlev, nlev,  nlev, nlevi, nlev, nlev,
                                 nlev, nlev, nlevi, nlev, nlev,  nlev, nlev};
  std::vector<const Real*> ptr_array = {wthv_sec, shoc_mix, u_wind,  v_wind, dz_zi, dz_zt, pres,
                                        brunt,    zt_grid,  zi_grid, tke,    tk,    tkh,   isotropy};

  // Sync to device
  ScreamDeepCopy::copy_to_device({obklen, pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  view_1d
    obklen_d(temp_1d_d[0]),
    pblh_d(temp_1d_d[1]);

  view_2d
    wthv_sec_d(temp_2d_d[0]),
    shoc_mix_d(temp_2d_d[1]),
    u_wind_d(temp_2d_d[2]),
    v_wind_d(temp_2d_d[3]),
    dz_zi_d(temp_2d_d[4]),
    dz_zt_d(temp_2d_d[5]),
    pres_d(temp_2d_d[6]),
    brunt_d(temp_2d_d[7]),
    zt_grid_d(temp_2d_d[8]),
    zi_grid_d(temp_2d_d[9]),
    tke_d(temp_2d_d[10]),
    tk_d(temp_2d_d[11]),
    tkh_d(temp_2d_d[12]),
    isotropy_d(temp_2d_d[13]);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);

  // Local variable workspace
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlevi_packs, 3, policy);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar obklen_s{obklen_d(i)};
    const Scalar pblh_s{pblh_d(i)};

    const auto wthv_sec_s = ekat::subview(wthv_sec_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);
    const auto u_wind_s = ekat::subview(u_wind_d, i);
    const auto v_wind_s = ekat::subview(v_wind_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto pres_s = ekat::subview(pres_d, i);
    const auto brunt_s = ekat::subview(brunt_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto tk_s = ekat::subview(tk_d, i);
    const auto tkh_s = ekat::subview(tkh_d, i);
    const auto isotropy_s = ekat::subview(isotropy_d, i);

    SHF::shoc_tke(team,nlev,nlevi,dtime,wthv_sec_s,shoc_mix_s,dz_zi_s,dz_zt_s,pres_s,
                  u_wind_s,v_wind_s,brunt_s,obklen_s,zt_grid_s,zi_grid_s,pblh_s,
                  workspace,
                  tke_s,tk_s,tkh_s,isotropy_s);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {tke_d, tk_d, tkh_d, isotropy_d};
  ekat::device_to_host({tke, tk, tkh, isotropy}, shcol, nlev, inout_views, true);
}

} // namespace shoc
} // namespace scream
