
! TODO: Add rules about ,value, reference scalars, and arrays in the interface


module cpp_interface_mod
  use params, only: crm_rknd, crm_iknd, crm_lknd
  use iso_c_binding
  implicit none

  public :: scream_session_init
  public :: scream_session_finalize

  interface

    subroutine crm(ncrms_in, pcols_in, dt_gl, plev, &
                   crm_input_bflxls, crm_input_wndls, &
                   crm_input_zmid, crm_input_zint, &
                   crm_input_pmid, crm_input_pint, crm_input_pdel, &
                   crm_input_ul, crm_input_vl, crm_input_tl, &
                   crm_input_qccl, crm_input_qiil, crm_input_ql, &
                   crm_input_tau00, crm_input_phis, crm_input_ps, &
#ifdef MMF_ESMT
                   crm_input_ul_esmt, crm_input_vl_esmt, &
#endif
                   crm_input_t_vt, crm_input_q_vt, &
                   crm_input_nccn, crm_input_nc_nuceat_tend, crm_input_ni_activated, &
                   crm_state_u_wind, crm_state_v_wind, crm_state_w_wind, &
                   crm_state_temperature, crm_state_qt, crm_state_qp, crm_state_qn, &
                   crm_state_qc, crm_state_nc, crm_state_qr, crm_state_nr, &
                   crm_state_qi, crm_state_ni, crm_state_qm, crm_state_bm, &
                   crm_state_t_prev, crm_state_q_prev, &
                   crm_rad_qrad, crm_rad_temperature, crm_rad_qv, crm_rad_qc, &
                   crm_rad_qi, crm_rad_cld, crm_rad_nc, crm_rad_ni, &
                   crm_output_subcycle_factor, &
                   crm_output_cld, crm_output_cldtop, crm_output_gicewp, crm_output_gliqwp, &
                   crm_output_mctot, crm_output_mcup, crm_output_mcdn, crm_output_mcuup, crm_output_mcudn, &
                   crm_output_qc_mean, crm_output_qi_mean, crm_output_qs_mean, crm_output_qg_mean, crm_output_qr_mean, &
                   crm_output_mu_crm, crm_output_md_crm, crm_output_eu_crm, crm_output_du_crm, crm_output_ed_crm, &
                   crm_output_flux_qt, crm_output_flux_u, crm_output_flux_v, crm_output_fluxsgs_qt, &
                   crm_output_tkez, crm_output_tkew, crm_output_tkesgsz, crm_output_tkz, &
                   crm_output_flux_qp, crm_output_precflux, crm_output_qt_trans, crm_output_qp_trans, &
                   crm_output_qp_fall, crm_output_qp_evp, crm_output_qp_src, crm_output_qt_ls, &
                   crm_output_t_ls, crm_output_jt_crm, crm_output_mx_crm, crm_output_cltot, &
                   crm_output_clhgh, crm_output_clmed, crm_output_cllow, &
                   crm_output_sltend, crm_output_qltend, crm_output_qcltend, crm_output_qiltend, &
                   crm_output_t_vt_tend, crm_output_q_vt_tend, crm_output_t_vt_ls, crm_output_q_vt_ls, &
#ifdef MMF_MOMENTUM_FEEDBACK
                   crm_output_ultend, crm_output_vltend, &
#endif 
                   crm_output_tk, crm_output_tkh, &
                   crm_output_qcl, crm_output_qci, crm_output_qpl, crm_output_qpi, &
                   crm_output_precc, crm_output_precl, crm_output_precsc, &
                   crm_output_precsl, crm_output_prec_crm,        & 
#ifdef MMF_ESMT
                   crm_output_u_tend_esmt, crm_output_v_tend_esmt, &
#endif
                   crm_clear_rh, &
                   lat0, long0, gcolp, igstep,  &
                   use_VT, VT_wn_max, &
                   microphysics_scheme, turbulence_scheme, &
                   use_crm_accel, crm_accel_factor, crm_accel_uv) bind(C,name="crm")
      use params, only: crm_rknd, crm_iknd, crm_lknd
      use iso_c_binding, only: c_bool, c_char
      implicit none
      logical(c_bool), value :: use_VT
      integer(crm_iknd), value :: VT_wn_max
      logical(c_bool), value :: use_crm_accel, crm_accel_uv
      integer(crm_iknd), value :: ncrms_in, pcols_in, plev, igstep
      real(crm_rknd), value :: dt_gl, crm_accel_factor
      integer(crm_iknd), dimension(*) :: gcolp
      character(kind=c_char) :: microphysics_scheme(*)
      character(kind=c_char) :: turbulence_scheme(*)

      real(crm_rknd), dimension(*) :: crm_input_bflxls, crm_input_wndls, &
                                      crm_input_zmid, crm_input_zint, &
                                      crm_input_pmid, crm_input_pint, crm_input_pdel, &
                                      crm_input_ul, crm_input_vl, crm_input_tl, &
                                      crm_input_qccl, crm_input_qiil, crm_input_ql, &
                                      crm_input_tau00, crm_input_phis, crm_input_ps, &
#ifdef MMF_ESMT
                                      crm_input_ul_esmt, crm_input_vl_esmt, &
#endif
                                      crm_input_t_vt, crm_input_q_vt, &
                                      crm_input_nccn, crm_input_nc_nuceat_tend, crm_input_ni_activated, &
                                      crm_state_u_wind, crm_state_v_wind, crm_state_w_wind, &
                                      crm_state_temperature, crm_state_qt, crm_state_qp, crm_state_qn, &
                                      crm_state_qc, crm_state_nc, crm_state_qr, crm_state_nr, &
                                      crm_state_qi, crm_state_ni, crm_state_qm, crm_state_bm, &
                                      crm_state_t_prev, crm_state_q_prev, &
                                      crm_rad_qrad, crm_rad_temperature, crm_rad_qv, crm_rad_qc, &
                                      crm_rad_qi, crm_rad_cld, crm_rad_nc, crm_rad_ni, &
                                      crm_output_subcycle_factor, &
                                      crm_output_cld, crm_output_cldtop, crm_output_gicewp, crm_output_gliqwp, &
                                      crm_output_mctot, crm_output_mcup, crm_output_mcdn, crm_output_mcuup, crm_output_mcudn, &
                                      crm_output_qc_mean, crm_output_qi_mean, crm_output_qs_mean, crm_output_qg_mean, crm_output_qr_mean, &
                                      crm_output_mu_crm, crm_output_md_crm, crm_output_eu_crm, crm_output_du_crm, crm_output_ed_crm, &
                                      crm_output_flux_qt, crm_output_flux_u, crm_output_flux_v, crm_output_fluxsgs_qt, &
                                      crm_output_tkez, crm_output_tkew, crm_output_tkesgsz, crm_output_tkz, &
                                      crm_output_flux_qp, crm_output_precflux, crm_output_qt_trans, crm_output_qp_trans, &
                                      crm_output_qp_fall, crm_output_qp_evp, crm_output_qp_src, &
                                      crm_output_qt_ls, crm_output_t_ls, &
                                      crm_output_jt_crm, crm_output_mx_crm, crm_output_cltot, &
                                      crm_output_clhgh, crm_output_clmed, crm_output_cllow, &
                                      crm_output_sltend, crm_output_qltend, crm_output_qcltend, crm_output_qiltend, &
                                      crm_output_t_vt_tend, crm_output_q_vt_tend, crm_output_t_vt_ls, crm_output_q_vt_ls, &
#ifdef MMF_MOMENTUM_FEEDBACK
                                      crm_output_ultend, crm_output_vltend, &
#endif
                                      crm_output_tk, crm_output_tkh, &
                                      crm_output_qcl, crm_output_qci, crm_output_qpl, crm_output_qpi, &
                                      crm_output_precc, crm_output_precl, crm_output_precsc, &
                                      crm_output_precsl, crm_output_prec_crm,        & 
#ifdef MMF_ESMT
                                      crm_output_u_tend_esmt, crm_output_v_tend_esmt, &
#endif
                                      crm_clear_rh, lat0, long0
    end subroutine crm


    subroutine setparm() bind(C,name="setparm")
      ! Do nothing
    end subroutine setparm

    subroutine scream_session_init() bind(C,name="scream_session_init")
      ! Do nothing
    end subroutine scream_session_init

    subroutine scream_session_finalize() bind(C,name="scream_session_finalize")
      ! Do nothing
    end subroutine scream_session_finalize

  end interface
end module cpp_interface_mod
