#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "share/scream_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "shoc_functions.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

// Glue functions to call fortran from from C++ with the Data struct
void shoc_grid                                      (SHOCGridData &d);
void shoc_diag_obklen                               (SHOCObklenData &d);
void update_host_dse                                (SHOCEnergydseData &d);
void shoc_energy_fixer                              (SHOCEnergyfixerData &d);
void shoc_energy_integrals                          (SHOCEnergyintData &d);
void shoc_energy_total_fixer                        (SHOCEnergytotData &d);
void shoc_energy_threshold_fixer                    (SHOCEnergythreshfixerData &d);
void shoc_energy_dse_fixer                          (SHOCEnergydsefixerData &d);
void calc_shoc_vertflux                             (SHOCVertfluxData &d);
void calc_shoc_varorcovar                           (SHOCVarorcovarData &d);
void compute_tmpi                                   (SHOCComptmpiData &d);
void dp_inverse                                     (SHOCDpinverseData &d);
void sfc_fluxes                                     (SHOCSfcfluxesData &d);
void impli_srf_stress_term                          (SHOCImplsrfstressData &d);
void tke_srf_flux_term                              (SHOCTkesrffluxData &d);
void integ_column_stability                         (SHOCColstabData &d);
void check_tke                                      (SHOCCheckTkeData &d);
void shoc_tke                                       (SHOCTkeData &d);
void compute_shr_prod                               (SHOCTkeshearData &d);
void isotropic_ts                                   (SHOCIsotropicData &d);
void adv_sgs_tke                                    (SHOCAdvsgstkeData &d);
void eddy_diffusivities                             (SHOCEddydiffData &d);
void shoc_length                                    (SHOCLengthData &d);
void compute_brunt_shoc_length                      (SHOCBruntlengthData &d);
void compute_l_inf_shoc_length                      (SHOCInflengthData &d);
void compute_conv_vel_shoc_length                   (SHOCConvvelData &d);
void compute_conv_time_shoc_length                  (SHOCConvtimeData &d);
void compute_shoc_mix_shoc_length                   (SHOCMixlengthData &d);
void check_length_scale_shoc_length                 (SHOCMixcheckData &d);
void fterms_input_for_diag_third_shoc_moment        (SHOCFterminputthirdmomsData &d);
void aa_terms_diag_third_shoc_moment                (SHOCAAdiagthirdmomsData &d);
void f0_to_f5_diag_third_shoc_moment                (SHOCFtermdiagthirdmomsData &d);
void omega_terms_diag_third_shoc_moment             (SHOCOmegadiagthirdmomsData &d);
void x_y_terms_diag_third_shoc_moment               (SHOCXYdiagthirdmomsData &d);
void w3_diag_third_shoc_moment                      (SHOCW3diagthirdmomsData &d);
void clipping_diag_third_shoc_moments               (SHOCClipthirdmomsData &d);
void shoc_diag_second_moments_srf                   (SHOCSecondMomentSrfData& d);
void linear_interp                                  (SHOCLinearInterpData &d);
void diag_third_shoc_moments                        (SHOCDiagThirdMomData &d);
void compute_diag_third_shoc_moment                 (SHOCCompThirdMomData &d);
void shoc_assumed_pdf                               (SHOCAssumedpdfData &d);
void shoc_assumed_pdf_tilde_to_real                 (SHOCPDFtildeData &d);
void shoc_assumed_pdf_vv_parameters                 (SHOCPDFvvparamData &d);
void shoc_assumed_pdf_thl_parameters                (SHOCPDFthlparamData &d);
void shoc_assumed_pdf_qw_parameters                 (SHOCPDFqwparamData &d);
void shoc_assumed_pdf_inplume_correlations          (SHOCPDFinplumeData &d);
void shoc_assumed_pdf_compute_temperature           (SHOCPDFcomptempData &d);
void shoc_assumed_pdf_compute_qs                    (SHOCPDFcompqsData &d);
void shoc_assumed_pdf_compute_s                     (SHOCPDFcompsData &d);
void shoc_assumed_pdf_compute_sgs_liquid            (SHOCPDFcompsgsliqData &d);
void shoc_assumed_pdf_compute_cloud_liquid_variance (SHOCPDFcompcloudvarData &d);
void shoc_assumed_pdf_compute_liquid_water_flux     (SHOCPDFcompliqfluxData &d);
void shoc_assumed_pdf_compute_buoyancy_flux         (SHOCPDFcompbuoyfluxData &d);
void shoc_diag_second_moments_ubycond               (SHOCSecondMomentUbycondData& d);
void shoc_pblintd_init_pot                          (SHOCPblintdInitPotData &d);
void shoc_pblintd_cldcheck                          (SHOCPblintdCldCheckData& d);
void diag_second_moments_lbycond                    (DiagSecondMomentsLbycondData& d);
void diag_second_moments                            (DiagSecondMomentsData& d);
void diag_second_shoc_moments                       (DiagSecondShocMomentsData& d);
void compute_shoc_vapor                             (ComputeShocVaporData& d);
void update_prognostics_implicit                    (UpdatePrognosticsImplicitData& d);

void shoc_main(ShocMainData& d);
void pblintd_height(PblintdHeightData& d);
extern "C" { // _f function decls

void calc_shoc_varorcovar_f(Int shcol, Int nlev, Int nlevi, Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
                            Real *invar1, Real *invar2, Real *varorcovar);
void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw,
                          Real* ustar2, Real* wstar);
void shoc_diag_second_moments_ubycond_f(Int shcol, Real* thl, Real* qw, Real* wthl,
                          Real* wqw, Real* qwthl, Real* uw, Real* vw, Real* wtke);
void update_host_dse_f(Int shcol, Int nlev, Real* thlm, Real* shoc_ql, Real* exner, Real* zt_grid,
                       Real* phis, Real* host_dse);
void compute_diag_third_shoc_moment_f(Int shcol, Int nlev, Int nlevi, Real* w_sec,
                                      Real* thl_sec, Real* wthl_sec, Real* tke,
                                      Real* dz_zt, Real* dz_zi, Real* isotropy_zi,
                                      Real* brunt_zi, Real* w_sec_zi, Real* thetal_zi,
                                      Real* w3);
void shoc_pblintd_init_pot_f(Int shcol, Int nlev, Real* thl, Real* ql, Real* q, Real* thv);
void compute_shoc_mix_shoc_length_f(Int nlev, Int shcol, Real* tke, Real* brunt,
                                    Real* tscale, Real* zt_grid, Real* l_inf, Real* shoc_mix);
void check_tke_f(Int shcol, Int nlev, Real* tke);
void linear_interp_f(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh);
void clipping_diag_third_shoc_moments_f(Int nlevi, Int shcol, Real *w_sec_zi,
                                        Real *w3);
void shoc_energy_integrals_f(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int);
void compute_brunt_shoc_length_f(Int nlev, Int nlevi, Int shcol, Real* dz_zt, Real* thv,
                                 Real* thv_zi, Real* brunt);
void compute_l_inf_shoc_length_f(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);
void check_length_scale_shoc_length_f(Int nlev, Int shcol, Real* host_dx, Real* host_dy,
                                      Real* shoc_mix);
void compute_conv_vel_shoc_length_f(Int nlev, Int shcol, Real *pblh, Real *zt_grid,
                                    Real *dz_zt, Real *thv, Real *wthv_sec,
                                    Real *conv_vel);
void diag_second_moments_lbycond_f(Int shcol, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar,
                                  Real* wthl_sec, Real* wqw_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* thl_sec,
                                  Real* qw_sec, Real* qwthl_sec);
void diag_second_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke, Real* isotropy,
                          Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix, Real* thl_sec, Real* qw_sec,
                          Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec);
void diag_second_shoc_moments_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke, Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec);
void shoc_diag_obklen_f(Int shcol, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc,
                        Real* thl_sfc, Real* cldliq_sfc, Real* qv_sfc, Real* ustar, Real* kbfs, Real* obklen);
void shoc_pblintd_cldcheck_f(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh);
void compute_conv_time_shoc_length_f(Int shcol, Real *pblh, Real *conv_vel,
                                     Real *tscale);
void shoc_length_f(Int shcol, Int nlev, Int nlevi, Real* host_dx, Real* host_dy, Real* pblh, Real* tke,
                   Real* zt_grid, Real* zi_grid, Real*dz_zt, Real* dz_zi, Real* wthv_sec, Real*thetal,
                   Real* thv, Real*brunt, Real* shoc_mix);
void shoc_energy_fixer_f(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* zt_grid,
                         Real* zi_grid, Real* se_b, Real* ke_b, Real* wv_b, Real* wl_b,
                         Real* se_a, Real* ke_a, Real* wv_a, Real* wl_a, Real* wthl_sfc,
                         Real* wqw_sfc, Real* rho_zt, Real* tke, Real* pint,
                         Real* host_dse);
void compute_shoc_vapor_f(Int shcol, Int nlev, Real* qw, Real* ql, Real* qv);
void update_prognostics_implicit_f(Int shcol, Int nlev, Int nlevi, Int num_tracer, Real dtime, Real* dz_zt, Real* dz_zi, Real* rho_zt, Real* zt_grid, Real* zi_grid, Real* tk, Real* tkh, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc, Real* wtracer_sfc, Real* thetal, Real* qw, Real* tracer, Real* tke, Real* u_wind, Real* v_wind);
void diag_third_shoc_moments_f(Int shcol, Int nlev, Int nlevi, Real* w_sec, Real* thl_sec,
                               Real* wthl_sec, Real* isotropy, Real* brunt, Real* thetal,
                               Real* tke, Real* dz_zt, Real* dz_zi, Real* zt_grid, Real* zi_grid,
                               Real* w3);
void shoc_assumed_pdf_f(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* w_field,
                        Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* w_sec, Real* wqw_sec,
                        Real* qwthl_sec, Real* w3, Real* pres, Real* zt_grid, Real* zi_grid,
                        Real* shoc_cldfrac, Real* shoc_ql, Real* wqls, Real* wthv_sec, Real* shoc_ql2);
void compute_tmpi_f(Int nlevi, Int shcol, Real dtime, Real *rho_zi, Real *dz_zi, Real *tmpi);
void dp_inverse_f(Int nlev, Int shcol, Real *rho_zt, Real *dz_zt, Real *rdp_zt);

void shoc_main_f(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* host_dx, Real* host_dy, Real* thv, Real* zt_grid, Real* zi_grid, Real* pres, Real* presi, Real* pdel, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* wtracer_sfc, Int num_qtracers, Real* w_field, Real* exner, Real* phis, Real* host_dse, Real* tke, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk, Real* shoc_ql, Real* shoc_cldfrac, Real* pblh, Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec, Real* qw_sec, Real* qwthl_sec, Real* wthl_sec, Real* wqw_sec, Real* wtke_sec, Real* uw_sec, Real* vw_sec, Real* w3, Real* wqls_sec, Real* brunt, Real* shoc_ql2);
void pblintd_height_f(Int shcol, Int nlev, Real* z, Real* u, Real* v, Real* ustar, Real* thv, Real* thv_ref, Real* pblh, Real* rino, bool* check);
} // end _f function decls

}  // namespace shoc
}  // namespace scream

#endif // SCREAM_SHOC_FUNCTIONS_F90_HPP
