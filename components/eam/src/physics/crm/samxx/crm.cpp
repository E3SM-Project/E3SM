
#include "pre_timeloop.h"
#include "post_timeloop.h"
#include "timeloop.h"
#include "vars.h"


extern "C" void crm(int ncrms_in, int pcols_in, real dt_gl, int plev, real *crm_input_bflxls_p, 
                    real *crm_input_wndls_p, real *crm_input_zmid_p, real *crm_input_zint_p, 
                    real *crm_input_pmid_p, real *crm_input_pint_p, real *crm_input_pdel_p, 
                    real *crm_input_ul_p, real *crm_input_vl_p, 
                    real *crm_input_tl_p, real *crm_input_qccl_p, real *crm_input_qiil_p, 
                    real *crm_input_ql_p, real *crm_input_tau00_p, 
                    real *crm_state_u_wind_p, real *crm_state_v_wind_p, real *crm_state_w_wind_p, 
                    real *crm_state_temperature_p, 
                    real *crm_state_qt_p, real *crm_state_qp_p, real *crm_state_qn_p, real *crm_rad_qrad_p, 
                    real *crm_rad_temperature_p, 
                    real *crm_rad_qv_p, real *crm_rad_qc_p, real *crm_rad_qi_p, real *crm_rad_cld_p, 
                    real *crm_output_subcycle_factor_p, 
                    real *crm_output_prectend_p, real *crm_output_precstend_p, real *crm_output_cld_p, 
                    real *crm_output_cldtop_p, 
                    real *crm_output_gicewp_p, real *crm_output_gliqwp_p, real *crm_output_mctot_p, 
                    real *crm_output_mcup_p, real *crm_output_mcdn_p, 
                    real *crm_output_mcuup_p, real *crm_output_mcudn_p, real *crm_output_qc_mean_p, 
                    real *crm_output_qi_mean_p, real *crm_output_qs_mean_p, 
                    real *crm_output_qg_mean_p, real *crm_output_qr_mean_p, real *crm_output_mu_crm_p, 
                    real *crm_output_md_crm_p, real *crm_output_eu_crm_p, 
                    real *crm_output_du_crm_p, real *crm_output_ed_crm_p, real *crm_output_flux_qt_p, 
                    real *crm_output_flux_u_p, real *crm_output_flux_v_p, 
                    real *crm_output_fluxsgs_qt_p, real *crm_output_tkez_p, real *crm_output_tkew_p, real *crm_output_tkesgsz_p, 
                    real *crm_output_tkz_p, real *crm_output_flux_qp_p, 
                    real *crm_output_precflux_p, real *crm_output_qt_trans_p, real *crm_output_qp_trans_p, 
                    real *crm_output_qp_fall_p, real *crm_output_qp_evp_p, 
                    real *crm_output_qp_src_p, real *crm_output_qt_ls_p, real *crm_output_t_ls_p, 
                    real *crm_output_jt_crm_p, real *crm_output_mx_crm_p, real *crm_output_cltot_p, 
                    real *crm_output_clhgh_p, real *crm_output_clmed_p, real *crm_output_cllow_p, 
                    real *crm_output_sltend_p, real *crm_output_qltend_p, real *crm_output_qcltend_p, real *crm_output_qiltend_p, 
#ifdef MMF_MOMENTUM_FEEDBACK
                    real *crm_output_ultend_p, real *crm_output_vltend_p,
#endif
                    real *crm_output_tk_p, real *crm_output_tkh_p, real *crm_output_qcl_p, 
                    real *crm_output_qci_p, real *crm_output_qpl_p, real *crm_output_qpi_p, 
                    real *crm_output_z0m_p, real *crm_output_taux_p, real *crm_output_tauy_p, real *crm_output_precc_p,
                    real *crm_output_precl_p, real *crm_output_precsc_p, 
                    real *crm_output_precsl_p, real *crm_output_prec_crm_p, real *crm_clear_rh_p,
                    real *lat0_p, real *long0_p, int *gcolp_p, 
                    int igstep_in,
                    bool use_crm_accel_in, real crm_accel_factor_in, bool crm_accel_uv_in) {

  dt_glob = dt_gl;
  pcols = pcols_in;
  ncrms = ncrms_in;
  igstep = igstep_in;
  use_crm_accel = use_crm_accel_in;
  crm_accel_factor = crm_accel_factor_in;
  crm_accel_uv = crm_accel_uv_in;

  create_and_copy_inputs(crm_input_bflxls_p, crm_input_wndls_p, crm_input_zmid_p, crm_input_zint_p, 
                         crm_input_pmid_p, crm_input_pint_p, crm_input_pdel_p, crm_input_ul_p, crm_input_vl_p, 
                         crm_input_tl_p, crm_input_qccl_p, crm_input_qiil_p, crm_input_ql_p, crm_input_tau00_p, 
                         crm_state_u_wind_p, crm_state_v_wind_p, crm_state_w_wind_p, crm_state_temperature_p, 
                         crm_state_qt_p, crm_state_qp_p, crm_state_qn_p, crm_rad_qrad_p, crm_output_subcycle_factor_p, 
                         lat0_p, long0_p, gcolp_p, crm_output_cltot_p, crm_output_clhgh_p, crm_output_clmed_p, 
                         crm_output_cllow_p);

  copy_outputs(crm_state_u_wind_p, crm_state_v_wind_p, crm_state_w_wind_p, crm_state_temperature_p, 
               crm_state_qt_p, crm_state_qp_p, crm_state_qn_p, crm_rad_temperature_p, 
               crm_rad_qv_p, crm_rad_qc_p, crm_rad_qi_p, crm_rad_cld_p, crm_output_subcycle_factor_p, 
               crm_output_prectend_p, crm_output_precstend_p, crm_output_cld_p, crm_output_cldtop_p, 
               crm_output_gicewp_p, crm_output_gliqwp_p, 
               crm_output_mctot_p, crm_output_mcup_p, crm_output_mcdn_p, 
               crm_output_mcuup_p, crm_output_mcudn_p, 
               crm_output_qc_mean_p, crm_output_qi_mean_p, crm_output_qs_mean_p, 
               crm_output_qg_mean_p, crm_output_qr_mean_p, crm_output_mu_crm_p, 
               crm_output_md_crm_p, crm_output_eu_crm_p, 
               crm_output_du_crm_p, crm_output_ed_crm_p, crm_output_flux_qt_p, 
               crm_output_flux_u_p, crm_output_flux_v_p, 
               crm_output_fluxsgs_qt_p, crm_output_tkez_p, crm_output_tkew_p, crm_output_tkesgsz_p, crm_output_tkz_p, 
               crm_output_flux_qp_p, crm_output_precflux_p, crm_output_qt_trans_p, crm_output_qp_trans_p, 
               crm_output_qp_fall_p, crm_output_qp_evp_p, crm_output_qp_src_p, crm_output_qt_ls_p, 
               crm_output_t_ls_p, crm_output_jt_crm_p, crm_output_mx_crm_p, 
               crm_output_cltot_p, crm_output_clhgh_p, crm_output_clmed_p, crm_output_cllow_p, 
               crm_output_sltend_p, crm_output_qltend_p, crm_output_qcltend_p, crm_output_qiltend_p, 
#ifdef MMF_MOMENTUM_FEEDBACK
               crm_output_ultend_p, crm_output_vltend_p,
#endif
               crm_output_tk_p, crm_output_tkh_p, 
               crm_output_qcl_p, crm_output_qci_p, crm_output_qpl_p, crm_output_qpi_p, 
               crm_output_z0m_p, crm_output_taux_p, crm_output_tauy_p, 
               crm_output_precc_p, crm_output_precl_p, crm_output_precsc_p, crm_output_precsl_p, 
               crm_output_prec_crm_p, crm_clear_rh_p);

  allocate();

  init_values();

  pre_timeloop();

  timeloop();

  post_timeloop();

  copy_outputs_and_destroy(crm_state_u_wind_p, crm_state_v_wind_p, crm_state_w_wind_p, crm_state_temperature_p, 
                           crm_state_qt_p, crm_state_qp_p, crm_state_qn_p, crm_rad_temperature_p, 
                           crm_rad_qv_p, crm_rad_qc_p, crm_rad_qi_p, crm_rad_cld_p, crm_output_subcycle_factor_p, 
                           crm_output_prectend_p, crm_output_precstend_p, crm_output_cld_p, crm_output_cldtop_p, 
                           crm_output_gicewp_p, crm_output_gliqwp_p, 
                           crm_output_mctot_p, crm_output_mcup_p, crm_output_mcdn_p, 
                           crm_output_mcuup_p, crm_output_mcudn_p, 
                           crm_output_qc_mean_p, crm_output_qi_mean_p, crm_output_qs_mean_p, 
                           crm_output_qg_mean_p, crm_output_qr_mean_p, crm_output_mu_crm_p, 
                           crm_output_md_crm_p, crm_output_eu_crm_p, 
                           crm_output_du_crm_p, crm_output_ed_crm_p, crm_output_flux_qt_p, 
                           crm_output_flux_u_p, crm_output_flux_v_p, 
                           crm_output_fluxsgs_qt_p, crm_output_tkez_p, crm_output_tkew_p, crm_output_tkesgsz_p, crm_output_tkz_p, 
                           crm_output_flux_qp_p, crm_output_precflux_p, crm_output_qt_trans_p, crm_output_qp_trans_p, 
                           crm_output_qp_fall_p, crm_output_qp_evp_p, crm_output_qp_src_p, crm_output_qt_ls_p, 
                           crm_output_t_ls_p, crm_output_jt_crm_p, crm_output_mx_crm_p, 
                           crm_output_cltot_p, crm_output_clhgh_p, crm_output_clmed_p, crm_output_cllow_p, 
                           crm_output_sltend_p, crm_output_qltend_p, crm_output_qcltend_p, crm_output_qiltend_p, 
#ifdef MMF_MOMENTUM_FEEDBACK
                           crm_output_ultend_p, crm_output_vltend_p,
#endif
                           crm_output_tk_p, crm_output_tkh_p, 
                           crm_output_qcl_p, crm_output_qci_p, crm_output_qpl_p, crm_output_qpi_p, 
                           crm_output_z0m_p, crm_output_taux_p, crm_output_tauy_p, 
                           crm_output_precc_p, crm_output_precl_p, crm_output_precsc_p, crm_output_precsl_p, 
                           crm_output_prec_crm_p, crm_clear_rh_p);

  finalize();
  
  yakl::fence();
}

