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

void update_host_dse_c(Int shcol, Int nlev, Real *thlm, Real *shoc_ql,
                       Real *exner, Real *zt_grid, Real *phis, Real *host_dse);

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

void shoc_length_c(Int shcol, Int nlev, Int nlevi, Real *tke, Real *host_dx,
                   Real *host_dy, Real *pblh, Real *zt_grid, Real *zi_grid,
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
void shoc_diag_second_moments_srf_c(Int shcol, Real* wthl, Real* uw, Real* vw,
                                   Real* ustar2, Real* wstar);

void diag_third_shoc_moments_c(Int shoc, Int nlev, Int nlevi, Real *w_sec, 
                               Real *thl_sec, Real *qw_sec, Real *qwthl_sec,
                               Real *wthl_sec, Real *isotropy, Real *brunt,
                               Real *thetal, Real *tke, Real *wthv_sec,
                               Real *dz_zt, Real *dz_zi, Real *zt_grid,
                               Real *zi_grid, Real *w3);

void compute_diag_third_shoc_moment_c(Int shcol, Int nlev, Int nlevi, Real *w_sec,
                                      Real *thl_sec, Real *qw_sec, Real *qwthl_sec,
                                      Real *wthl_sec, Real *tke, Real *dz_zt, 
                                      Real *dz_zi, Real *zt_grid, Real *zi_grid,
                                      Real *isotropy_zi, Real *brunt_zi, 
                                      Real *w_sec_zi, Real *thetal_zi, 
                                      Real *wthv_sec_zi, Real *w3);

void linear_interp_c(Real *x1, Real *x2, Real *y1, Real *y2, Int km1,
                     Int km2, Int ncol, Real minthresh);
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

void shoc_diag_second_moments_ubycond_c(Int shcol, Real* thl, Real* qw, Real* wthl,
                                       Real* wqw, Real* qwthl, Real* uw, Real* vw,
                                       Real* wtke);

void shoc_pblintd_init_pot_c(Int shcol, Int nlev, Real* thl, Real* ql, Real* q, Real* thv);

}

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

void shoc_grid(SHOCGridData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol(), d.nlev(), d.nlevi(), d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void update_host_dse(SHOCEnergydseData &d) {
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  update_host_dse_c(d.shcol(), d.nlev(), d.thlm, d.shoc_ql, d.exner,
                    d.zt_grid, d.phis, d.host_dse);
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
  shoc_length_c(d.shcol(),d.nlev(),d.nlevi(),d.tke,d.host_dx,d.host_dy,
                d.pblh,d.zt_grid,d.zi_grid,d.dz_zt,d.dz_zi,d.thetal,
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
  shoc_diag_second_moments_srf_c(d.shcol(), d.wthl, d.uw, d.vw, d.ustar2, d.wstar);
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
                                   d.qw_sec,d.qwthl_sec,d.wthl_sec,d.tke,d.dz_zt,
                                   d.dz_zi,d.zt_grid,d.zi_grid,d.isotropy_zi,
                                   d.brunt_zi,d.w_sec_zi,d.thetal_zi,d.wthv_sec_zi,
                                   d.w3);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void linear_interp(SHOCLinearintData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::TransposeDirection::c2f>();
  linear_interp_c(d.x1,d.x2,d.y1,d.y2,d.nlev(),d.nlevi(),d.shcol(),d.minthresh);
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
  shoc_diag_second_moments_ubycond_c(d.shcol(), d.thl, d.qw, d.wthl, d.wqw, d.qwthl, d.uw, d.vw, d.wtke);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void shoc_pblintd_init_pot(SHOCPblintdInitPotData& d)
{
  shoc_init(d.nlev(), true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_pblintd_init_pot_c(d.shcol(), d.nlev(), d.thl, d.ql, d.q, d.thv);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

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

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw, Real* ustar2, Real* wstar)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Pack1      = typename ekat::Pack<Real, 1>;
  using view_1d    = typename SHOC::view_1d<Pack1>;

  Kokkos::Array<view_1d, 3> temp_d;
  ekat::host_to_device({wthl, uw, vw}, shcol, temp_d);

  // inputs
  view_1d
    wthl_d (temp_d[0]),
    uw_d   (temp_d[1]),
    vw_d   (temp_d[2]);

  // outputs
  view_1d ustar2_d("ustar2", shcol),
          wstar_d ("wstar", shcol);

  Kokkos::parallel_for("parallel_moments_srf", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar wthl_s{wthl_d(i)[0]};
     Scalar uw_s{uw_d(i)[0]};
     Scalar vw_s{vw_d(i)[0]};

     Scalar ustar2_s{0};
     Scalar wstar_s{0};

     SHOC::shoc_diag_second_moments_srf(wthl_s, uw_s, vw_s, ustar2_s, wstar_s);

     ustar2_d(i)[0] = ustar2_s;
     wstar_d(i)[0]  = wstar_s;
   });

  Kokkos::Array<view_1d, 2> out_views = {ustar2_d, wstar_d};
  ekat::device_to_host({ustar2, wstar}, shcol, out_views);
}

void shoc_diag_second_moments_ubycond_f(Int shcol, Real* thl, Real* qw, Real* wthl, Real* wqw, Real* qwthl, Real* uw, Real* vw,
      Real* wtke)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Pack1      = typename ekat::Pack<Real, 1>;
  using view_1d    = typename SHOC::view_1d<Pack1>;

  view_1d thl_d  ("thl"  ,shcol),
          qw_d   ("qw"   ,shcol),
          qwthl_d("qwthl",shcol),
          wthl_d ("wthl" ,shcol),
          wqw_d  ("wqw"  ,shcol),
          uw_d   ("uw"   ,shcol),
          vw_d   ("vw"   ,shcol),
          wtke_d ("wtke" ,shcol);

  Kokkos::parallel_for("parallel_moments_ubycond", shcol, KOKKOS_LAMBDA (const int& i) {

    Scalar thl_s{0.};
    Scalar qw_s{0.};
    Scalar wthl_s{0.};
    Scalar wqw_s{0.};
    Scalar qwthl_s{0.};
    Scalar uw_s{0.};
    Scalar vw_s{0.};
    Scalar wtke_s{0.};

    SHOC::shoc_diag_second_moments_ubycond(thl_s, qw_s, wthl_s, wqw_s, qwthl_s, uw_s, vw_s, wtke_s);

    thl_d(i)[0]   = thl_s;
    qw_d(i)[0]    = qw_s;
    wthl_d(i)[0]  = wthl_s;
    wqw_d(i)[0]   = wqw_s;
    qwthl_d(i)[0] = qwthl_s;
    uw_d(i)[0]    = uw_s;
    vw_d(i)[0]    = vw_s;
    wtke_d(i)[0]  = wtke_s;

  });

  Kokkos::Array<view_1d, 8> host_views = {thl_d, qw_d, qwthl_d, wthl_d, wqw_d, uw_d, vw_d, wtke_d};

  ekat::device_to_host({thl, qw, qwthl, wthl, wqw, uw, vw, wtke}, shcol, host_views);
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
  ekat::pack::host_to_device({thl, ql, q}, shcol, nlev, temp_d, true);

  view_2d thl_d(temp_d[0]),
          ql_d (temp_d[1]),
          q_d  (temp_d[2]);

  view_2d thv_d("thv", shcol, nlev);

  const Int nlev_pack = ekat::pack::npack<Spack>(nlev);
  const auto policy = ekat::util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto thl_1d = ekat::util::subview(thl_d, i);
    const auto ql_1d  = ekat::util::subview(ql_d, i);
    const auto q_1d   = ekat::util::subview(q_d, i);
    const auto thv_1d = ekat::util::subview(thv_d, i);

    SHOC::shoc_pblintd_init_pot(team, nlev, thl_1d, ql_1d, q_1d, thv_1d);
  });

  Kokkos::Array<view_2d, 1> inout_views = {thv_d};
  ekat::pack::device_to_host({thv}, {shcol}, {nlev}, inout_views, true);
}

} // namespace shoc
} // namespace scream
