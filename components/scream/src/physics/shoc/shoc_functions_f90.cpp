#include "shoc_functions_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "shoc_f90.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to SHOC fortran calls. The stubs below will link to fortran definitions in shoc_iso_c.f90
//

extern "C" {

void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                 Real zvir, Real latvap, Real latice, Real karman,
                 Real* pref_mid, int nbot_shoc, int ntop_shoc);

void shoc_grid_c(int shcol, int nlev, int nlevi, Real *zt_grid, Real *zi_grid,
                 Real *pdel, Real *dz_zt, Real *dzi_zi, Real *rho_zt);

void shoc_diag_obklen_c(Int shcol, Real *uw_sfc, Real *vw_sfc, Real *wthl_sfc,
                        Real *wqw_sfc, Real *thl_sfc, Real *cldliq_sfc,
                        Real *qv_sfc, Real *ustar, Real *kbfs, Real *obklen);

void update_host_dse_c(Int shcol, Int nlev, Real *thlm, Real *shoc_ql,
                       Real *exner, Real *zt_grid, Real *phis, Real *host_dse);

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

void sfc_fluxes_c(Int shcol, Real dtime, Real *rho_zi_sfc, Real *rdp_zt_sfc,
                  Real *wthl_sfc, Real *wqw_sfc, Real *wtke_sfc, Real *thetal,
                  Real *qw, Real *tke);

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
                   Real *host_dy, Real *pblh, Real *tke, Real *zt_grid, Real *zi_grid,
                   Real *dz_zt, Real *dz_zi, Real *thetal, Real *wthv_sec,
                   Real *thv, Real *brunt, Real *shoc_mix);

void compute_brunt_shoc_length_c(Int nlev, Int nlevi, Int shcol ,Real *dz_zt,
                                 Real *thv, Real *thv_zi, Real *brunt);

void compute_l_inf_shoc_length_c(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);

void compute_conv_vel_shoc_length_c(Int nlev, Int shcol, Real *pblh, Real *zt_grid,
                                    Real *dz_zt, Real *thv, Real *wthv_sec,
                                    Real *conv_vel);

void compute_conv_time_shoc_length_c(Int shcol, Real *pblh, Real *conv_vel,
                                     Real *tscale);

void compute_shoc_mix_shoc_length_c(Int nlev, Int shcol, Real *tke, Real* brunt,
                                    Real *tscale, Real *zt_grid, Real *l_inf,
				    Real *shoc_mix);

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
				    Real wthl_sec_kb, Real thl_sec, Real thl_sec_kc,
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
                               Real *thl_sec, Real *qw_sec, Real *qwthl_sec,
                               Real *wthl_sec, Real *isotropy, Real *brunt,
                               Real *thetal, Real *tke, Real *wthv_sec,
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

void shoc_assumed_pdf_tilda_to_real_c(Real w_first, Real sqrtw2, Real* w1);

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

} // end _c function decls

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

void calc_shoc_varorcovar(SHOCVarorcovarData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol(), d.nlev(), d.nlevi(), d.tunefac, d.isotropy_zi, d.tkh_zi,
                         d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_tmpi(SHOCComptmpiData &d){
  shoc_init(d.nlevi()-1, true); // nlev=nlevi-1
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_tmpi_c(d.nlevi(), d.shcol(), d.dtime, d.rho_zi, d.dz_zi, d.tmpi);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void dp_inverse(SHOCDpinverseData &d){
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  dp_inverse_c(d.nlev(), d.shcol(), d.rho_zt, d.dz_zt, d.rdp_zt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void sfc_fluxes(SHOCSfcfluxesData &d){
  shoc_init(1, true); // single layer function
  d.transpose<ekat::TransposeDirection::c2f>();
  sfc_fluxes_c(d.shcol(), d.dtime, d.rho_zi_sfc, d.rdp_zt_sfc, d.wthl_sfc,
               d.wqw_sfc, d.wtke_sfc, d.thetal, d.qw, d.tke);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void impli_srf_stress_term(SHOCImplsrfstressData &d){
  shoc_init(1, true); // single layer function
  d.transpose<ekat::TransposeDirection::c2f>();
  impli_srf_stress_term_c(d.shcol(), d.rho_zi_sfc, d.uw_sfc, d.vw_sfc,
                          d.u_wind_sfc, d.v_wind_sfc, d.ksrf);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void tke_srf_flux_term(SHOCTkesrffluxData &d){
  shoc_init(1, true); // single layer function
  d.transpose<ekat::TransposeDirection::c2f>();
  tke_srf_flux_term_c(d.shcol(), d.uw_sfc, d.vw_sfc, d.wtke_sfc);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_grid(SHOCGridData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol(), d.nlev(), d.nlevi(), d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_diag_obklen(SHOCObklenData &d){
  shoc_init(1, true); // single level function
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_diag_obklen_c(d.shcol(), d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc,
                     d.thl_sfc, d.cldliq_sfc, d.qv_sfc, d.ustar, d.kbfs,
                     d.obklen);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void update_host_dse(SHOCEnergydseData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  update_host_dse_c(d.shcol(), d.nlev(), d.thlm, d.shoc_ql, d.exner,
                    d.zt_grid, d.phis, d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_fixer(SHOCEnergyfixerData &d){
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_fixer_c(d.shcol(), d.nlev(), d.nlevi(), d.dtime, d.nadv,
                      d.zt_grid, d.zi_grid, d.se_b, d.ke_b, d.wv_b,
                      d.wl_b, d.se_a, d.ke_a, d.wv_a, d.wl_a, d.wthl_sfc,
                      d.wqw_sfc, d.rho_zt, d.tke, d.pint,
                      d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_integrals(SHOCEnergyintData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_integrals_c(d.shcol(), d.nlev(), d.host_dse, d.pdel,
                          d.rtm, d.rcm, d.u_wind, d.v_wind,
                          d.se_int, d.ke_int, d.wv_int, d.wl_int);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_total_fixer(SHOCEnergytotData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_total_fixer_c(d.shcol(), d.nlev(), d.nlevi(), d.dtime, d.nadv,
                            d.zt_grid, d.zi_grid,
                            d.se_b, d.ke_b, d.wv_b, d.wl_b,
                            d.se_a, d.ke_a, d.wv_a, d.wl_a,
                            d.wthl_sfc, d.wqw_sfc, d.rho_zt,
                            d.te_a, d.te_b);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_threshold_fixer(SHOCEnergythreshfixerData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_threshold_fixer_c(d.shcol(), d.nlev(), d.nlevi(),
                          d.pint, d.tke, d.te_a, d.te_b,
			  d.se_dis, d.shoctop);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_energy_dse_fixer(SHOCEnergydsefixerData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_energy_dse_fixer_c(d.shcol(), d.nlev(),
                          d.se_dis, d.shoctop, d.host_dse);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void calc_shoc_vertflux(SHOCVertfluxData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol(), d.nlev(), d.nlevi(), d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void check_tke(SHOCCheckTkeData &d){
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  check_tke_c(d.shcol(), d.nlev(), d.tke);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_tke(SHOCTkeData &d){
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_tke_c(d.shcol(), d.nlev(), d.nlevi(), d.dtime, d.wthv_sec, d.shoc_mix,
              d.dz_zi, d.dz_zt, d.pres, d.u_wind, d.v_wind, d.brunt, d.obklen,
	      d.zt_grid, d.zi_grid, d.pblh, d.tke, d.tk, d.tkh, d.isotropy);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void integ_column_stability(SHOCColstabData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  integ_column_stability_c(d.nlev(), d.shcol(), d.dz_zt, d.pres, d.brunt, d.brunt_int);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_shr_prod(SHOCTkeshearData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_shr_prod_c(d.nlevi(), d.nlev(), d.shcol(), d.dz_zi, d.u_wind,
                       d.v_wind, d.sterm);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void isotropic_ts(SHOCIsotropicData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  isotropic_ts_c(d.nlev(), d.shcol(), d.brunt_int, d.tke, d.a_diss,
                 d.brunt, d.isotropy);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void adv_sgs_tke(SHOCAdvsgstkeData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  adv_sgs_tke_c(d.nlev(), d.shcol(), d.dtime, d.shoc_mix, d.wthv_sec,
                d.sterm_zt, d.tk, d.tke, d.a_diss);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void eddy_diffusivities(SHOCEddydiffData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  eddy_diffusivities_c(d.nlev(), d.shcol(), d.obklen, d.pblh, d.zt_grid,
     d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_length(SHOCLengthData &d){
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_length_c(d.shcol(),d.nlev(),d.nlevi(),d.host_dx,d.host_dy,
                d.pblh,d.tke,d.zt_grid,d.zi_grid,d.dz_zt,d.dz_zi,d.thetal,
                d.wthv_sec,d.thv,d.brunt,d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_brunt_shoc_length(SHOCBruntlengthData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_brunt_shoc_length_c(d.nlev(),d.nlevi(),d.shcol(),d.dz_zt,d.thv,d.thv_zi,d.brunt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_l_inf_shoc_length(SHOCInflengthData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_l_inf_shoc_length_c(d.nlev(),d.shcol(),d.zt_grid,d.dz_zt,d.tke,d.l_inf);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_conv_vel_shoc_length(SHOCConvvelData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_conv_vel_shoc_length_c(d.nlev(),d.shcol(),d.pblh,d.zt_grid,
                                 d.dz_zt,d.thv,d.wthv_sec,d.conv_vel);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_conv_time_shoc_length(SHOCConvtimeData &d) {
  shoc_init(42, true); // fake nlev
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_conv_time_shoc_length_c(d.shcol(),d.pblh,d.conv_vel,d.tscale);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_shoc_mix_shoc_length(SHOCMixlengthData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_shoc_mix_shoc_length_c(d.nlev(),d.shcol(),d.tke,d.brunt,d.tscale,
                                 d.zt_grid,d.l_inf,d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void check_length_scale_shoc_length(SHOCMixcheckData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  check_length_scale_shoc_length_c(d.nlev(),d.shcol(),d.host_dx,d.host_dy,d.shoc_mix);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void clipping_diag_third_shoc_moments(SHOCClipthirdmomsData &d) {
  shoc_init(d.nlevi()-1, true); // nlev = nlevi - 1
  d.transpose<ekat::TransposeDirection::c2f>();
  clipping_diag_third_shoc_moments_c(d.nlevi(),d.shcol(),d.w_sec_zi,d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void aa_terms_diag_third_shoc_moment(SHOCAAdiagthirdmomsData &d){
  shoc_init(1, true);
  aa_terms_diag_third_shoc_moment_c(d.omega0,d.omega1,d.omega2,d.x0,d.x1,d.y0,d.y1,
                                    &d.aa0,&d.aa1);
}

void fterms_input_for_diag_third_shoc_moment(SHOCFterminputthirdmomsData &d){
  shoc_init(1, true);
  fterms_input_for_diag_third_shoc_moment_c(d.dz_zi,d.dz_zt,d.dz_zt_kc,
                                     d.isotropy_zi,d.brunt_zi,d.thetal_zi,
				     &d.thedz,&d.thedz2,&d.iso,
				     &d.isosqrd,&d.buoy_sgs2,&d.bet2);
}

void f0_to_f5_diag_third_shoc_moment(SHOCFtermdiagthirdmomsData &d){
  shoc_init(1, true);
  f0_to_f5_diag_third_shoc_moment_c(d.thedz,d.thedz2,d.bet2,d.iso,d.isosqrd,
                                    d.wthl_sec,d.wthl_sec_kc,d.wthl_sec_kb,
                                    d.thl_sec,d.thl_sec_kc,d.thl_sec_kb,
                                    d.w_sec,d.w_sec_kc,d.w_sec_zi,
                                    d.tke,d.tke_kc,
                                    &d.f0,&d.f1,&d.f2,&d.f3,&d.f4,&d.f5);
}

void omega_terms_diag_third_shoc_moment(SHOCOmegadiagthirdmomsData &d){
  shoc_init(1, true);
  omega_terms_diag_third_shoc_moment_c(d.buoy_sgs2,d.f3,d.f4,
                                       &d.omega0,&d.omega1,&d.omega2);
}

void x_y_terms_diag_third_shoc_moment(SHOCXYdiagthirdmomsData &d){
  shoc_init(1, true);
  x_y_terms_diag_third_shoc_moment_c(d.buoy_sgs2,d.f0,d.f1,d.f2,
                                     &d.x0,&d.y0,&d.x1,&d.y1);
}

void w3_diag_third_shoc_moment(SHOCW3diagthirdmomsData &d){
  shoc_init(1, true);
  w3_diag_third_shoc_moment_c(d.aa0,d.aa1,d.x0,d.x1,d.f5,&d.w3);
}

void shoc_diag_second_moments_srf(SHOCSecondMomentSrfData& d)
{
  shoc_init(42, true); // fake nlev
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_diag_second_moments_srf_c(d.shcol(), d.wthl_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_third_shoc_moments(SHOCDiagThirdMomData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_third_shoc_moments_c(d.shcol(),d.nlev(),d.nlevi(),d.w_sec,d.thl_sec,d.qw_sec,
                            d.qwthl_sec,d.wthl_sec,d.isotropy,d.brunt,d.thetal,
                            d.tke,d.wthv_sec,d.dz_zt,d.dz_zi,d.zt_grid,d.zi_grid,
                            d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void compute_diag_third_shoc_moment(SHOCCompThirdMomData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  compute_diag_third_shoc_moment_c(d.shcol(),d.nlev(),d.nlevi(),d.w_sec,d.thl_sec,
                                   d.wthl_sec,d.tke,d.dz_zt,d.dz_zi,d.isotropy_zi,
                                   d.brunt_zi,d.w_sec_zi,d.thetal_zi,d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void linear_interp(SHOCLinearInterpData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  linear_interp_c(d.x1,d.x2,d.y1,d.y2,d.nlev(),d.nlevi(),d.shcol(),d.minthresh);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_assumed_pdf(SHOCAssumedpdfData &d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_assumed_pdf_c(d.shcol(), d.nlev(), d.nlevi(), d.thetal, d.qw, d.w_field,
                     d.thl_sec, d.qw_sec, d.wthl_sec, d.w_sec, d.wqw_sec,
                     d.qwthl_sec, d.w3, d.pres, d.zt_grid, d.zi_grid,
                     d.shoc_cldfrac, d.shoc_ql, d.wqls, d.wthv_sec, d.shoc_ql2);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_assumed_pdf_tilda_to_real(SHOCPDFtildaData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_tilda_to_real_c(d.w_first, d.sqrtw2, &d.w1);
}

void shoc_assumed_pdf_vv_parameters(SHOCPDFvvparamData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_vv_parameters_c(d.w_first,d.w_sec,d.w3var,
                                   &d.Skew_w,&d.w1_1,&d.w1_2,&d.w2_1,&d.w2_2,&d.a);
}

void shoc_assumed_pdf_thl_parameters(SHOCPDFthlparamData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_thl_parameters_c(d.wthlsec,d.sqrtw2,d.sqrtthl,d.thlsec,d.thl_first,
                                    d.w1_1,d.w1_2,d.Skew_w,d.a,d.dothetal_skew,
                                    &d.thl1_1,&d.thl1_2,&d.thl2_1,&d.thl2_2,&d.sqrtthl2_1,
                                    &d.sqrtthl2_2);
}

void shoc_assumed_pdf_qw_parameters(SHOCPDFqwparamData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_qw_parameters_c(d.wqwsec,d.sqrtw2,d.Skew_w,d.sqrtqt,d.qwsec,
                                    d.w1_1,d.w1_2,d.qw_first,d.a,
                                    &d.qw1_1,&d.qw1_2,&d.qw2_1,&d.qw2_2,&d.sqrtqw2_1,
                                    &d.sqrtqw2_2);
}

void shoc_assumed_pdf_inplume_correlations(SHOCPDFinplumeData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_inplume_correlations_c(d.sqrtqw2_1,d.sqrtthl2_1,d.a,
                                          d.sqrtqw2_2,d.sqrtthl2_2,
                                          d.qwthlsec,d.qw1_1,d.qw_first,d.thl1_1,
			                  d.thl_first,d.qw1_2,d.thl1_2,
                                          &d.r_qwthl_1);
}

void shoc_assumed_pdf_compute_temperature(SHOCPDFcomptempData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_temperature_c(d.thl1, d.basepres, d.pval, &d.Tl1);
}

void shoc_assumed_pdf_compute_qs(SHOCPDFcompqsData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_qs_c(d.Tl1_1,d.Tl1_2,d.pval,
                                &d.qs1,&d.beta1,&d.qs2,&d.beta2);
}

void shoc_assumed_pdf_compute_s(SHOCPDFcompsData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_s_c(d.qw1,d.qs1,d.beta,d.pval,d.thl2,d.qw2,
                               d.sqrtthl2,d.sqrtqw2,d.r_qwthl,
			       &d.s,&d.std_s,&d.qn,&d.C);
}

void shoc_assumed_pdf_compute_sgs_liquid(SHOCPDFcompsgsliqData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_sgs_liquid_c(d.a, d.ql1, d.ql2, &d.shoc_ql);
}

void shoc_assumed_pdf_compute_cloud_liquid_variance(SHOCPDFcompcloudvarData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_cloud_liquid_variance_c(d.a,d.s1,d.ql1,d.C1,
                                   d.std_s1,d.s2,d.ql2,d.C2,d.std_s2,d.shoc_ql,
                                   &d.shoc_ql2);
}

void shoc_assumed_pdf_compute_liquid_water_flux(SHOCPDFcompliqfluxData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_liquid_water_flux_c(d.a,d.w1_1,d.w_first,d.ql1,
                                               d.w1_2,d.ql2,&d.wqls);
}

void shoc_assumed_pdf_compute_buoyancy_flux(SHOCPDFcompbuoyfluxData &d)
{
  shoc_init(1, true);
  shoc_assumed_pdf_compute_buoyancy_flux_c(d.wthlsec,d.epsterm,d.wqwsec,
                                           d.pval,d.wqls,&d.wthv_sec);
}

void shoc_diag_second_moments_ubycond(SHOCSecondMomentUbycondData& d)
{
  shoc_init(42, true); // Fake nlev
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_diag_second_moments_ubycond_c(d.shcol(), d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_pblintd_init_pot(SHOCPblintdInitPotData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_pblintd_init_pot_c(d.shcol(), d.nlev(), d.thl, d.ql, d.q, d.thv);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_second_moments_lbycond(DiagSecondMomentsLbycondData& d)
{
  shoc_init(1, true);  // Single level routine
  diag_second_moments_lbycond_c(d.shcol(), d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar, d.wthl_sec, d.wqw_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.thl_sec, d.qw_sec, d.qwthl_sec);
}

void diag_second_moments(DiagSecondMomentsData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_second_moments_c(d.shcol(), d.nlev(), d.nlevi(), d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh, d.tk, d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.w_sec);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void diag_second_shoc_moments(DiagSecondShocMomentsData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  diag_second_shoc_moments_c(d.shcol(), d.nlev(), d.nlevi(), d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy, d.tkh, d.tk, d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.thl_sec, d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.w_sec);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_pblintd_cldcheck(SHOCPblintdCldCheckData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_pblintd_cldcheck_c(d.shcol(), d.nlev(), d.nlevi(), d.zi, d.cldn, d.pblh);
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

  Kokkos::Array<view_2d, num_arrays> temp_d;
  Kokkos::Array<int, num_arrays> dim1_sizes     = {shcol,       shcol,  shcol, shcol,  shcol,  shcol};
  Kokkos::Array<int, num_arrays> dim2_sizes     = {nlevi,       nlevi,  nlevi, nlev,   nlev,   nlevi};
  Kokkos::Array<const Real*, num_arrays> ptr_array = {isotropy_zi, tkh_zi, dz_zi, invar1, invar2, varorcovar};

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
  Kokkos::Array<view_2d, 1> inout_views = {varorcovar_d};
  ekat::device_to_host<int,1>({varorcovar}, {shcol}, {nlevi}, inout_views, true);
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

  Kokkos::Array<view_2d, num_arrays> temp_d;
  Kokkos::Array<Int, num_arrays> dim1_sizes     = {shcol,  shcol, shcol, shcol};
  Kokkos::Array<Int, num_arrays> dim2_sizes     = {nlevi,  nlevi, nlev,  nlevi};
  Kokkos::Array<const Real*, num_arrays> ptr_array = {tkh_zi, dz_zi, invar, vertflux};

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
  Kokkos::Array<view_2d, 1> inout_views = {vertflux_d};
  ekat::device_to_host<Int,1>({vertflux}, {{shcol}}, {{nlevi}}, inout_views, true);
}

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Pack1      = typename ekat::Pack<Real, 1>;
  using view_1d    = typename SHOC::view_1d<Pack1>;

  Kokkos::Array<view_1d, 3> temp_d;
  ekat::host_to_device({wthl_sfc, uw_sfc, vw_sfc}, shcol, temp_d);

  // inputs
  view_1d
    wthl_sfc_d (temp_d[0]),
    uw_sfc_d   (temp_d[1]),
    vw_sfc_d   (temp_d[2]);

  // outputs
  view_1d ustar2_d("ustar2", shcol),
          wstar_d ("wstar", shcol);

  Kokkos::parallel_for("parallel_moments_srf", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar wthl_sfc_s{wthl_sfc_d(i)[0]};
     Scalar uw_sfc_s{uw_sfc_d(i)[0]};
     Scalar vw_sfc_s{vw_sfc_d(i)[0]};

     Scalar ustar2_s{0};
     Scalar wstar_s{0};

     SHOC::shoc_diag_second_moments_srf(wthl_sfc_s, uw_sfc_s, vw_sfc_s, ustar2_s, wstar_s);

     ustar2_d(i)[0] = ustar2_s;
     wstar_d(i)[0]  = wstar_s;
   });

  Kokkos::Array<view_1d, 2> out_views = {ustar2_d, wstar_d};
  ekat::device_to_host({ustar2, wstar}, shcol, out_views);
}

void shoc_diag_second_moments_ubycond_f(Int shcol, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec,
      Real* wtke_sec)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Pack1      = typename ekat::Pack<Real, 1>;
  using view_1d    = typename SHOC::view_1d<Pack1>;

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

    thl_sec_d(i)[0]   = thl_sec_s;
    qw_sec_d(i)[0]    = qw_sec_s;
    wthl_sec_d(i)[0]  = wthl_sec_s;
    wqw_sec_d(i)[0]   = wqw_sec_s;
    qwthl_sec_d(i)[0] = qwthl_sec_s;
    uw_sec_d(i)[0]    = uw_sec_s;
    vw_sec_d(i)[0]    = vw_sec_s;
    wtke_sec_d(i)[0]  = wtke_sec_s;

  });

  Kokkos::Array<view_1d, 8> host_views = {thl_sec_d, qw_sec_d, qwthl_sec_d, wthl_sec_d, wqw_sec_d, uw_sec_d, vw_sec_d, wtke_sec_d};

  ekat::device_to_host({thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec}, shcol, host_views);
}

void update_host_dse_f(Int shcol, Int nlev, Real* thlm, Real* shoc_ql, Real* exner, Real* zt_grid,
                       Real* phis, Real* host_dse)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 1> temp_1d_d;
  Kokkos::Array<view_2d, 5> temp_2d_d;
  Kokkos::Array<int, 5> dim1_sizes     = {shcol,  shcol, shcol, shcol, shcol};
  Kokkos::Array<int, 5> dim2_sizes     = {nlev,  nlev, nlev,  nlev, nlev};
  Kokkos::Array<const Real*, 5> ptr_array = {thlm, shoc_ql, exner, zt_grid, host_dse};

  // Sync to device
  ekat::host_to_device({phis}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  view_1d phis_d(temp_1d_d[0]);

  view_2d
    thlm_d    (temp_2d_d[0]),
    shoc_ql_d (temp_2d_d[1]),
    exner_d   (temp_2d_d[2]),
    zt_grid_d (temp_2d_d[3]),
    host_dse_d(temp_2d_d[4]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar phis_s{phis_d(i)[0]};
    const auto thlm_s   = ekat::subview(thlm_d, i);
    const auto shoc_ql_s    = ekat::subview(shoc_ql_d, i);
    const auto exner_s    = ekat::subview(exner_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto host_dse_s = ekat::subview(host_dse_d, i);

    SHF::update_host_dse(team, nlev, thlm_s, shoc_ql_s, exner_s, zt_grid_s, phis_s, host_dse_s);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {host_dse_d};
  ekat::device_to_host<int,1>({host_dse}, {shcol}, {nlev}, inout_views, true);
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

  Kokkos::Array<view_2d, 11> temp_d;
  Kokkos::Array<size_t, 11> dim1_sizes     = {shcol,       shcol,
                                              shcol,       shcol,
                                              shcol,       shcol,
                                              shcol,       shcol,
                                              shcol,       shcol,
                                              shcol};
  Kokkos::Array<size_t, 11> dim2_sizes     = {nlev,        nlevi,
                                              nlevi,       nlev,
                                              nlev,        nlevi,
                                              nlevi,       nlevi,
                                              nlevi,       nlevi,
                                              nlevi};
  Kokkos::Array<const Real*, 11> ptr_array = {w_sec,       thl_sec,
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
  Kokkos::Array<view_2d, 1> inout_views = {w3_d};
  ekat::device_to_host<int,1>({w3}, {shcol}, {nlevi}, inout_views, true);
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

  Kokkos::Array<view_2d, num_arrays> temp_d;
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

  Kokkos::Array<view_2d, 1> inout_views = {thv_d};
  ekat::device_to_host<int,1>({thv}, {shcol}, {nlev}, inout_views, true);
}

void compute_shoc_mix_shoc_length_f(Int nlev, Int shcol, Real* tke, Real* brunt, Real* tscale,
                                    Real* zt_grid, Real* l_inf, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

   using Scalar     = typename SHF::Scalar;
   using Spack      = typename SHF::Spack;
   using Pack1d     = typename ekat::Pack<Real,1>;
   using view_1d    = typename SHF::view_1d<Pack1d>;
   using view_2d    = typename SHF::view_2d<Spack>;
   using KT         = typename SHF::KT;
   using ExeSpace   = typename KT::ExeSpace;
   using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 2> temp_1d_d;
  Kokkos::Array<view_2d, 4> temp_2d_d;
  Kokkos::Array<size_t, 4> dim1_sizes     = {shcol, shcol, shcol,   shcol};
  Kokkos::Array<size_t, 4> dim2_sizes     = {nlev,  nlev,  nlev,    nlev};
  Kokkos::Array<const Real*, 4> ptr_array = {tke,   brunt, zt_grid, shoc_mix};

  // Sync to device
  ekat::host_to_device({tscale,l_inf}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  view_1d
    tscale_d (temp_1d_d[0]),
    l_inf_d  (temp_1d_d[1]);

  view_2d
    tke_d     (temp_2d_d[0]),
    brunt_d   (temp_2d_d[1]),
    zt_grid_d (temp_2d_d[2]),
    shoc_mix_d  (temp_2d_d[3]);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar tscale_s{tscale_d(i)[0]};
    const Scalar l_inf_s {l_inf_d(i)[0]};

    const auto tke_s      = ekat::subview(tke_d, i);
    const auto brunt_s    = ekat::subview(brunt_d, i);
    const auto zt_grid_s  = ekat::subview(zt_grid_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    SHF::compute_shoc_mix_shoc_length(team, nlev, tke_s, brunt_s, tscale_s, zt_grid_s, l_inf_s,
                                      shoc_mix_s);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {shoc_mix_d};
  ekat::device_to_host<int,1>({shoc_mix}, {shcol}, {nlev}, inout_views, true);
}

void check_tke_f(Int shcol, Int nlev, Real* tke)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  Kokkos::Array<view_2d, 1> temp_2d_d;

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
  Kokkos::Array<view_2d, 1> inout_views = {tke_d};
  ekat::device_to_host<int,1>({tke}, {shcol}, {nlev}, inout_views, true);
}

void linear_interp_f(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_2d, 3> temp_2d_d;
  Kokkos::Array<size_t, 3> dim1_sizes     = {ncol, ncol, ncol};
  Kokkos::Array<size_t, 3> dim2_sizes     = {km1,  km2,  km1};
  Kokkos::Array<const Real*, 3> ptr_array = {x1,   x2,   y1};

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
  Kokkos::Array<view_2d, 1> inout_views = {y2_d};
  ekat::device_to_host<int,1>({y2}, {ncol}, {km2}, inout_views, true);
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

  Kokkos::Array<view_2d, 2> temp_d;
  Kokkos::Array<size_t, 2> dim1_sizes     = {shcol, shcol};
  Kokkos::Array<size_t, 2> dim2_sizes     = {nlevi, nlevi};
  Kokkos::Array<const Real*, 2> ptr_array = {w_sec_zi, w3};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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
  Kokkos::Array<view_2d, 1> inout_views = {w3_d};
  ekat::device_to_host<int,1>({w3}, {shcol}, {nlevi}, inout_views, true);
}

void shoc_energy_integrals_f(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_2d, 6> temp_d;
  Kokkos::Array<int, 6> dim1_sizes        = {shcol,   shcol, shcol, shcol, shcol,  shcol};
  Kokkos::Array<int, 6> dim2_sizes        = {nlev,     nlev, nlev,  nlev,  nlev,   nlev};
  Kokkos::Array<const Real*, 6> ptr_array = {host_dse, pdel, rtm,   rcm,   u_wind, v_wind};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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

    se_int_d(i)[0] = se_int_s;
    ke_int_d(i)[0] = ke_int_s;
    wv_int_d(i)[0] = wv_int_s;
    wl_int_d(i)[0] = wl_int_s;
  });

  // Sync back to host
  Kokkos::Array<view_1d, 4> inout_views = {se_int_d, ke_int_d, wv_int_d, wl_int_d};
  ekat::device_to_host<int,4>({se_int,ke_int,wv_int,wl_int},shcol,inout_views);
}

void diag_second_moments_lbycond_f(Int shcol, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar, Real* wthl_sec, Real* wqw_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* thl_sec, Real* qw_sec, Real* qwthl_sec)
{
  // TODO
}
void diag_second_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke, Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec)
{
  // TODO
}
void diag_second_shoc_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke, Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec)
{
  // TODO
}

void compute_brunt_shoc_length_f(Int nlev, Int nlevi, Int shcol, Real* dz_zt, Real* thv, Real* thv_zi, Real* brunt)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_2d, 4> temp_d;
  Kokkos::Array<int, 4> dim1_sizes        = {shcol, shcol, shcol,  shcol};
  Kokkos::Array<int, 4> dim2_sizes        = {nlev,  nlev,  nlevi,  nlev};
  Kokkos::Array<const Real*, 4> ptr_array = {dz_zt, thv,   thv_zi, brunt};

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
  Kokkos::Array<view_2d, 1> inout_views = {brunt_d};
  ekat::device_to_host<int,1>({brunt}, {shcol}, {nlev}, inout_views, true);
}

void compute_l_inf_shoc_length_f(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_2d, 3> temp_d;
  Kokkos::Array<int, 3> dim1_sizes        = {shcol,   shcol, shcol};
  Kokkos::Array<int, 3> dim2_sizes        = {nlev,     nlev, nlev};
  Kokkos::Array<const Real*, 3> ptr_array = {zt_grid, dz_zt, tke};

  // Sync to device
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

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

    l_inf_d(i)[0] = l_inf_s;
  });

  // Sync back to host
  Kokkos::Array<view_1d, 1> inout_views = {l_inf_d};
  ekat::device_to_host<int,1>({l_inf},shcol,inout_views);
}

void check_length_scale_shoc_length_f(Int nlev, Int shcol, Real* host_dx, Real* host_dy, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 2> temp_1d_d;
  Kokkos::Array<view_2d, 1> temp_2d_d;
  Kokkos::Array<int, 1> dim1_sizes        = {shcol};
  Kokkos::Array<int, 1> dim2_sizes        = {nlev};
  Kokkos::Array<const Real*, 1> ptr_array = {shoc_mix};

  // Sync to device
  ekat::host_to_device({host_dx,host_dy}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  view_1d
    host_dx_d(temp_1d_d[0]),
    host_dy_d(temp_1d_d[1]);

  view_2d
    shoc_mix_d(temp_2d_d[0]);
  
  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar host_dx_s{host_dx_d(i)[0]};
    const Scalar host_dy_s{host_dy_d(i)[0]};
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    SHF::check_length_scale_shoc_length(team, nlev, host_dx_s, host_dy_s, shoc_mix_s);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {shoc_mix_d};
  ekat::device_to_host<int,1>({shoc_mix}, {shcol}, {nlev}, inout_views, true);
}

void compute_conv_vel_shoc_length_f(Int nlev, Int shcol, Real *pblh, Real *zt_grid,
                                    Real *dz_zt, Real *thv, Real *wthv_sec,
                                    Real *conv_vel)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 1> temp_1d_d;
  Kokkos::Array<view_2d, 4> temp_2d_d;
  Kokkos::Array<int, 4> dim1_sizes        = {shcol,   shcol, shcol, shcol};
  Kokkos::Array<int, 4> dim2_sizes        = {nlev,     nlev, nlev,  nlev};
  Kokkos::Array<const Real*, 4> ptr_array = {zt_grid, dz_zt, thv,   wthv_sec};

  // Sync to device
  ekat::host_to_device({pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  // inputs
  view_1d
    pblh_d (temp_1d_d[0]);

  view_2d
    zt_grid_d (temp_2d_d[0]),
    dz_zt_d   (temp_2d_d[1]),
    thv_d     (temp_2d_d[2]),
    wthv_sec_d(temp_2d_d[3]);

  // outputs
  view_1d
    conv_vel_d("conv_vel", shcol);

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Inputs
    const Scalar pblh_s{pblh_d(i)[0]};
    const auto zt_grid_s  = ekat::subview(zt_grid_d, i);
    const auto dz_zt_s    = ekat::subview(dz_zt_d, i);
    const auto thv_s      = ekat::subview(thv_d, i);
    const auto wthv_sec_s = ekat::subview(wthv_sec_d, i);

    // Output
    Scalar conv_vel_s{0};

    SHF::compute_conv_vel_shoc_length(team, nlev, pblh_s, zt_grid_s, dz_zt_s, thv_s, wthv_sec_s,
                                      conv_vel_s);

    conv_vel_d(i)[0] = conv_vel_s;
  });

  // Sync back to host
  Kokkos::Array<view_1d, 1> inout_views = {conv_vel_d};
  ekat::device_to_host<int,1>({conv_vel},shcol,inout_views);
}

void shoc_diag_obklen_f(Int shcol, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc, Real* thl_sfc,
                        Real* cldliq_sfc, Real* qv_sfc, Real* ustar, Real* kbfs, Real* obklen)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 7> temp_d;
  Kokkos::Array<const Real*, 7> ptr_array = {uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, thl_sfc,
                                             cldliq_sfc, qv_sfc};

  // Sync to device
  ekat::host_to_device(ptr_array, shcol, temp_d);

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
    Scalar uw_sfc_s{uw_sfc_d(i)[0]};
    Scalar vw_sfc_s{vw_sfc_d(i)[0]};
    Scalar wthl_sfc_s{wthl_sfc_d(i)[0]};
    Scalar wqw_sfc_s{wqw_sfc_d(i)[0]};
    Scalar thl_sfc_s{thl_sfc_d(i)[0]};
    Scalar cldliq_sfc_s{cldliq_sfc_d(i)[0]};
    Scalar qv_sfc_s{qv_sfc_d(i)[0]};

    Scalar ustar_s{0};
    Scalar kbfs_s{0};
    Scalar obklen_s{0};

    SHF::shoc_diag_obklen(uw_sfc_s, vw_sfc_s, wthl_sfc_s, wqw_sfc_s, thl_sfc_s, cldliq_sfc_s, qv_sfc_s,
                          ustar_s, kbfs_s, obklen_s);

    ustar_d(i)[0] = ustar_s;
    kbfs_d(i)[0] = kbfs_s;
    obklen_d(i)[0] = obklen_s;
  });

  // Sync back to host
  Kokkos::Array<view_1d, 3> inout_views = {ustar_d, kbfs_d, obklen_d};
  ekat::device_to_host<int,3>({ustar, kbfs, obklen}, shcol, inout_views);
}

void shoc_pblintd_cldcheck_f(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh) {
  using SHOC    = Functions<Real, DefaultDevice>;
  using Pack1   = typename ekat::Pack<Real, 1>;
  using Scalar  = typename SHOC::Scalar;
  using view_2d = typename SHOC::view_2d<Pack1>;
  using view_1d = typename SHOC::view_1d<Pack1>;

  Kokkos::Array<size_t, 2> dim1  = {shcol, shcol};
  Kokkos::Array<size_t, 2> dim2  = {nlevi,  nlev};

  Kokkos::Array<view_2d, 2> cldcheck_2d;
  ekat::host_to_device({zi, cldn}, dim1, dim2, cldcheck_2d, true);

  view_2d
         zi_2d  (cldcheck_2d[0]),
         cldn_2d(cldcheck_2d[1]);

  Kokkos::Array<view_1d, 1> cldcheck_1d;
  ekat::host_to_device({pblh}, shcol, cldcheck_1d);

  view_1d pblh_1d (cldcheck_1d[0]);

  Kokkos::parallel_for("pblintd_cldcheck", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar zi_s   = zi_2d(i, nlev-1)[0];
     Scalar cldn_s = cldn_2d(i, nlev-1)[0];
     Scalar pblh_s = pblh_1d(i)[0];

     SHOC::shoc_pblintd_cldcheck(zi_s, cldn_s, pblh_s);

     pblh_1d(i)[0] = pblh_s;

  });

  Kokkos::Array<view_1d, 1> host_views = {pblh_1d};

  ekat::device_to_host<int,1>({pblh}, shcol, host_views);
}

void compute_conv_time_shoc_length_f(Int shcol, Real *pblh, Real *conv_vel, Real *tscale)
{
  using SHF       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHF::Scalar;
  using Pack1      = typename ekat::Pack<Real, 1>;
  using view_1d    = typename SHF::view_1d<Pack1>;

  Kokkos::Array<view_1d, 3> temp_d;
  ekat::host_to_device({pblh, conv_vel, tscale}, shcol, temp_d);

  view_1d
    pblh_d(temp_d[0]),
    conv_vel_d(temp_d[1]),
    tscale_d(temp_d[2]);

  Kokkos::parallel_for("compute_conv_time_shoc_length", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar pblh_s{pblh_d(i)[0]};
     Scalar conv_vel_s{conv_vel_d(i)[0]};
     Scalar tscale_s{tscale_d(i)[0]};

     SHF::compute_conv_time_shoc_length(pblh_s, conv_vel_s, tscale_s);

     conv_vel_d(i)[0] = conv_vel_s;
     tscale_d(i)[0]  = tscale_s;
   });

  Kokkos::Array<view_1d, 2> inout_views = {conv_vel_d, tscale_d};
  ekat::device_to_host({conv_vel, tscale}, shcol, inout_views);
}

void shoc_length_f(Int shcol, Int nlev, Int nlevi, Real* host_dx, Real* host_dy, Real* pblh, Real* tke,
                   Real* zt_grid, Real* zi_grid, Real*dz_zt, Real* dz_zi, Real* wthv_sec, Real*thetal,
                   Real* thv, Real*brunt, Real* shoc_mix)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Scalar     = typename SHF::Scalar;
  using Spack      = typename SHF::Spack;
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 3> temp_1d_d;
  Kokkos::Array<view_2d, 10> temp_2d_d;
  Kokkos::Array<int, 10> dim1_sizes = {shcol, shcol, shcol, shcol, shcol,
                                       shcol, shcol, shcol, shcol, shcol};
  Kokkos::Array<int, 10> dim2_sizes = {nlev, nlev, nlevi, nlev, nlevi,
                                       nlev, nlev, nlev,  nlev, nlev};
  Kokkos::Array<const Real*, 10> ptr_array = {tke,      zt_grid, zi_grid, dz_zt, dz_zi,
                                              wthv_sec, thetal,  thv,     brunt, shoc_mix};
  // Sync to device
  ekat::host_to_device({host_dx, host_dy, pblh}, shcol, temp_1d_d);
  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_2d_d, true);

  // inputs
  view_1d
    host_dx_d(temp_1d_d[0]),
    host_dy_d(temp_1d_d[1]),
    pblh_d(temp_1d_d[2]);

  view_2d
    tke_d(temp_2d_d[0]),
    zt_grid_d(temp_2d_d[1]),
    zi_grid_d(temp_2d_d[2]),
    dz_zt_d(temp_2d_d[3]),
    dz_zi_d(temp_2d_d[4]),
    wthv_sec_d(temp_2d_d[5]),
    thetal_d(temp_2d_d[6]),
    thv_d(temp_2d_d[7]),
    brunt_d(temp_2d_d[8]),
    shoc_mix_d(temp_2d_d[9]);

  // Local variable
  view_2d thv_zi_d("thv_zi", shcol, ekat::npack<Spack>(nlevi));

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Inputs
    const Scalar host_dx_s{host_dx_d(i)[0]};
    const Scalar host_dy_s{host_dy_d(i)[0]};
    const Scalar pblh_s{pblh_d(i)[0]};

    const auto tke_s = ekat::subview(tke_d, i);
    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto dz_zt_s = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s = ekat::subview(dz_zi_d, i);
    const auto wthv_sec_s = ekat::subview(wthv_sec_d, i);
    const auto thetal_s = ekat::subview(thetal_d, i);
    const auto thv_s = ekat::subview(thv_d, i);
    const auto thv_zi_s = ekat::subview(thv_zi_d, i);
    const auto brunt_s = ekat::subview(brunt_d, i);
    const auto shoc_mix_s = ekat::subview(shoc_mix_d, i);

    SHF::shoc_length(team,nlev,nlevi,host_dx_s,host_dy_s,pblh_s,tke_s,
                     zt_grid_s,zi_grid_s,dz_zt_s,dz_zi_s,wthv_sec_s,
                     thetal_s,thv_s,thv_zi_s,brunt_s,shoc_mix_s);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 2> out_views = {brunt_d,shoc_mix_d};
  ekat::device_to_host<int,2>({brunt,shoc_mix},shcol,nlev,out_views,true);
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
  using Pack1d     = typename ekat::Pack<Real,1>;
  using view_1d    = typename SHF::view_1d<Pack1d>;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  Kokkos::Array<view_1d, 10> temp_1d_d;
  Kokkos::Array<const Real*, 10> ptr_array_1d = {se_b, ke_b, wv_b, wl_b,     se_a,
                                                 ke_a, wv_a, wl_a, wthl_sfc, wqw_sfc};
  Kokkos::Array<view_2d, 6> temp_2d_d;
  Kokkos::Array<int, 6> dim1_sizes           = {shcol,   shcol,   shcol,
                                                shcol,   shcol,   shcol};
  Kokkos::Array<int, 6> dim2_sizes           = {nlev,    nlevi,   nlevi,
                                                nlev,    nlev,    nlev};
  Kokkos::Array<const Real*, 6> ptr_array_2d = {zt_grid, zi_grid, pint,
                                                rho_zt,  tke,     host_dse};

  // Sync to device
  ekat::host_to_device(ptr_array_1d, shcol, temp_1d_d);
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

  // Local variable
  view_2d rho_zi_d("rho_zi", shcol, ekat::npack<Spack>(nlevi));

  const Int nk_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar se_b_s{se_b_d(i)[0]};
    const Scalar ke_b_s{ke_b_d(i)[0]};
    const Scalar wv_b_s{wv_b_d(i)[0]};
    const Scalar wl_b_s{wl_b_d(i)[0]};
    const Scalar se_a_s{se_a_d(i)[0]};
    const Scalar ke_a_s{ke_a_d(i)[0]};
    const Scalar wv_a_s{wv_a_d(i)[0]};
    const Scalar wl_a_s{wl_a_d(i)[0]};
    const Scalar wthl_sfc_s{wthl_sfc_d(i)[0]};
    const Scalar wqw_sfc_s{wqw_sfc_d(i)[0]};

    const auto zt_grid_s = ekat::subview(zt_grid_d, i);
    const auto zi_grid_s = ekat::subview(zi_grid_d, i);
    const auto pint_s = ekat::subview(pint_d, i);
    const auto rho_zt_s = ekat::subview(rho_zt_d, i);
    const auto tke_s = ekat::subview(tke_d, i);
    const auto rho_zi_s = ekat::subview(rho_zi_d, i);
    const auto host_dse_s = ekat::subview(host_dse_d, i);

    SHF::shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid_s,zi_grid_s,se_b_s,
                           ke_b_s,wv_b_s,wl_b_s,se_a_s,ke_a_s,wv_a_s,wl_a_s,
                           wthl_sfc_s,wqw_sfc_s,rho_zt_s,tke_s,pint_s,rho_zi_s,host_dse_s);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {host_dse_d};
  ekat::device_to_host<int,1>({host_dse}, {shcol}, {nlev}, inout_views, true);
}

} // namespace shoc
} // namespace scream
