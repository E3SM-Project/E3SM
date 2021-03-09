
#pragma once

#include "samxx_const.h"


void allocate();


void init_values();


void finalize();


inline void perturb(real1d &arr, double mag) {
  for (int i=0; i<arr.get_totElems(); i++) {
    double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    arr.data()[i] *= (1.0 + r*mag);
  }
}


inline real arraySum(real1d &arr) {
  real sum=0;
  for (int i=0; i<arr.get_totElems(); i++) {
    sum += arr.data()[i];
  }
  return sum;
}


//////////////////////////////////////////////////////////////////////////////////
// These arrays use non-1 lower bounds in the Fortran code
// They must be indexed differently in the C++ code
//////////////////////////////////////////////////////////////////////////////////
extern real4d u             ; // Index as u             (    k , offy_u    +j , offx_u    +i , icrm )
extern real4d v             ; // Index as v             (    k , offy_v    +j , offx_v    +i , icrm )
extern real4d w             ; // Index as w             (    k , offy_w    +j , offx_w    +i , icrm )
extern real4d t             ; // Index as t             (    k , offy_s    +j , offx_s    +i , icrm )
extern real4d p             ; // Index as p             (    k , offy_p    +j , offx_p    +i , icrm )
extern real4d tke2          ; // Index as tke2          (    k , offy_s +j , offx_s +i , icrm )
extern real4d tk2           ; // Index as tk2           (    k , offy_tk2(or offy_d)  +j , offx_d  +i , icrm )
extern real3d sstxy         ; // Index as sstxy         (        offy_sstxy+j , offx_sstxy+i , icrm )
extern real2d fcory         ; // Index as fcory         (        offy_fcory+j                , icrm )
extern real5d sgs_field     ; // Index as sgs_field     (l , k , offy_s    +j , offx_s    +i , icrm )
extern real5d sgs_field_diag; // Index as sgs_field_diag(l , k , offy_d    +j , offx_d    +i , icrm )
extern real5d micro_field   ; // Index as sgs_field_diag(l , k , offy_s    +j , offx_s    +i , icrm )


void perturb_arrays();


void create_and_copy_inputs(real *crm_input_bflxls_p, real *crm_input_wndls_p, real *crm_input_zmid_p, real *crm_input_zint_p, 
                            real *crm_input_pmid_p, real *crm_input_pint_p, real *crm_input_pdel_p, real *crm_input_ul_p, real *crm_input_vl_p, 
                            real *crm_input_tl_p, real *crm_input_qccl_p, real *crm_input_qiil_p, real *crm_input_ql_p, real *crm_input_tau00_p, 
                            real *crm_state_u_wind_p, real *crm_state_v_wind_p, real *crm_state_w_wind_p, real *crm_state_temperature_p, 
                            real *crm_state_qt_p, real *crm_state_qp_p, real *crm_state_qn_p, real *crm_rad_qrad_p, real *crm_output_subcycle_factor_p, 
                            real *lat0_p, real *long0_p, int *gcolp_p, real *crm_output_cltot_p, real *crm_output_clhgh_p, real *crm_output_clmed_p,
                            real *crm_output_cllow_p);
                            


void copy_outputs(real *crm_state_u_wind_p, real *crm_state_v_wind_p, real *crm_state_w_wind_p, real *crm_state_temperature_p, 
                  real *crm_state_qt_p, real *crm_state_qp_p, real *crm_state_qn_p, real *crm_rad_temperature_p, 
                  real *crm_rad_qv_p, real *crm_rad_qc_p, real *crm_rad_qi_p, real *crm_rad_cld_p, real *crm_output_subcycle_factor_p, 
                  real *crm_output_prectend_p, real *crm_output_precstend_p, real *crm_output_cld_p, real *crm_output_cldtop_p, 
                  real *crm_output_gicewp_p, real *crm_output_gliqwp_p, real *crm_output_mctot_p, real *crm_output_mcup_p, real *crm_output_mcdn_p, 
                  real *crm_output_mcuup_p, real *crm_output_mcudn_p, real *crm_output_qc_mean_p, real *crm_output_qi_mean_p, real *crm_output_qs_mean_p, 
                  real *crm_output_qg_mean_p, real *crm_output_qr_mean_p, real *crm_output_mu_crm_p, real *crm_output_md_crm_p, real *crm_output_eu_crm_p, 
                  real *crm_output_du_crm_p, real *crm_output_ed_crm_p, real *crm_output_flux_qt_p, real *crm_output_flux_u_p, real *crm_output_flux_v_p, 
                  real *crm_output_fluxsgs_qt_p, real *crm_output_tkez_p, real *crm_output_tkew_p, real *crm_output_tkesgsz_p, real *crm_output_tkz_p, real *crm_output_flux_qp_p, 
                  real *crm_output_precflux_p, real *crm_output_qt_trans_p, real *crm_output_qp_trans_p, real *crm_output_qp_fall_p, real *crm_output_qp_evp_p, 
                  real *crm_output_qp_src_p, real *crm_output_qt_ls_p, real *crm_output_t_ls_p, real *crm_output_jt_crm_p, real *crm_output_mx_crm_p, real *crm_output_cltot_p, 
                  real *crm_output_clhgh_p, real *crm_output_clmed_p, real *crm_output_cllow_p, 
                  real *crm_output_sltend_p, real *crm_output_qltend_p, real *crm_output_qcltend_p, real *crm_output_qiltend_p, 
#ifdef MMF_MOMENTUM_FEEDBACK
                  real *crm_output_ultend_p, real *crm_output_vltend_p,
#endif
                  real *crm_output_tk_p, real *crm_output_tkh_p, real *crm_output_qcl_p, real *crm_output_qci_p, real *crm_output_qpl_p, real *crm_output_qpi_p, 
                  real *crm_output_z0m_p, real *crm_output_taux_p, real *crm_output_tauy_p, real *crm_output_precc_p, real *crm_output_precl_p, real *crm_output_precsc_p, 
                  real *crm_output_precsl_p, real *crm_output_prec_crm_p, real *crm_clear_rh_p);
                            


void copy_outputs_and_destroy(real *crm_state_u_wind_p, real *crm_state_v_wind_p, real *crm_state_w_wind_p, real *crm_state_temperature_p, 
                              real *crm_state_qt_p, real *crm_state_qp_p, real *crm_state_qn_p, real *crm_rad_temperature_p, 
                              real *crm_rad_qv_p, real *crm_rad_qc_p, real *crm_rad_qi_p, real *crm_rad_cld_p, real *crm_output_subcycle_factor_p, 
                              real *crm_output_prectend_p, real *crm_output_precstend_p, real *crm_output_cld_p, real *crm_output_cldtop_p, 
                              real *crm_output_gicewp_p, real *crm_output_gliqwp_p, real *crm_output_mctot_p, real *crm_output_mcup_p, real *crm_output_mcdn_p, 
                              real *crm_output_mcuup_p, real *crm_output_mcudn_p, real *crm_output_qc_mean_p, real *crm_output_qi_mean_p, real *crm_output_qs_mean_p, 
                              real *crm_output_qg_mean_p, real *crm_output_qr_mean_p, real *crm_output_mu_crm_p, real *crm_output_md_crm_p, real *crm_output_eu_crm_p, 
                              real *crm_output_du_crm_p, real *crm_output_ed_crm_p, real *crm_output_flux_qt_p, real *crm_output_flux_u_p, real *crm_output_flux_v_p, 
                              real *crm_output_fluxsgs_qt_p, real *crm_output_tkez_p, real *crm_output_tkew_p, real *crm_output_tkesgsz_p, real *crm_output_tkz_p, real *crm_output_flux_qp_p, 
                              real *crm_output_precflux_p, real *crm_output_qt_trans_p, real *crm_output_qp_trans_p, real *crm_output_qp_fall_p, real *crm_output_qp_evp_p, 
                              real *crm_output_qp_src_p, real *crm_output_qt_ls_p, real *crm_output_t_ls_p, real *crm_output_jt_crm_p, real *crm_output_mx_crm_p, real *crm_output_cltot_p, 
                              real *crm_output_clhgh_p, real *crm_output_clmed_p, real *crm_output_cllow_p, 
                              real *crm_output_sltend_p, real *crm_output_qltend_p, real *crm_output_qcltend_p, real *crm_output_qiltend_p, 
#ifdef MMF_MOMENTUM_FEEDBACK
                              real *crm_output_ultend_p, real *crm_output_vltend_p, 
#endif
                              real *crm_output_tk_p, real *crm_output_tkh_p, real *crm_output_qcl_p, real *crm_output_qci_p, real *crm_output_qpl_p, real *crm_output_qpi_p, 
                              real *crm_output_z0m_p, real *crm_output_taux_p, real *crm_output_tauy_p, real *crm_output_precc_p, real *crm_output_precl_p, real *crm_output_precsc_p, 
                              real *crm_output_precsl_p, real *crm_output_prec_crm_p, real *crm_clear_rh_p);


                            
extern int pcols;
extern int ncrms;

extern int  nstep                    ;
extern int  ncycle                   ;
extern int  icycle                   ;
extern int  na, nb, nc               ;
extern real at, bt, ct               ;
extern real dtn                      ;
extern real dtfactor                 ;
extern int  rank                     ;
extern int  ranknn                   ;
extern int  rankss                   ;
extern int  rankee                   ;
extern int  rankww                   ;
extern int  rankne                   ;
extern int  ranknw                   ;
extern int  rankse                   ;
extern int  ranksw                   ;
extern bool dompi                    ;
extern bool masterproc               ;
extern bool dostatis                 ;
extern bool dostatisrad              ;
extern int  nstatis                  ;
extern bool compute_reffc            ;
extern bool compute_reffi            ;
extern bool notopened2D              ;
extern bool notopened3D              ;
extern bool notopenedmom             ;
extern real dx                       ;
extern real dy                       ;
extern bool doconstdz                ;
extern int  nstop                    ;
extern int  nelapse                  ;
extern real dt                       ;
extern real day0                     ;
extern int  nrad                     ;
extern int  nrestart                 ;
extern int  nstat                    ;
extern int  nstatfrq                 ;
extern bool restart_sep              ;
extern int  nrestart_skip            ;
extern bool output_sep               ;
extern bool doisccp                  ;
extern bool domodis                  ;
extern bool domisr                   ;
extern bool dosimfilesout            ;
extern bool doSAMconditionals        ;
extern bool dosatupdnconditionals    ;
extern bool doscamiopdata            ;
extern bool dozero_out_day0          ;
extern int  nsave3Dstart             ;
extern int  nsave3Dend               ;
extern bool save3Dbin                ;
extern bool save3Dsep                ;
extern real qnsave3D                 ;
extern bool dogzip3D                 ;
extern bool rad3Dout                 ;
extern int  nsave2D                  ;
extern int  nsave2Dstart             ;
extern int  nsave2Dend               ;
extern bool save2Dbin                ;
extern bool save2Dsep                ;
extern bool save2Davg                ;
extern bool dogzip2D                 ;
extern int  nstatmom                 ;
extern int  nstatmomstart            ;
extern int  nstatmomend              ;
extern bool savemomsep               ;
extern bool savemombin               ;
extern int  nmovie                   ;
extern int  nmoviestart              ;
extern int  nmovieend                ;
extern bool isInitialized_scamiopdata;
extern bool wgls_holds_omega         ;

extern bool dosubsidence    ;
extern real ug              ;
extern real vg              ;
extern bool les             ;
extern bool sfc_flx_fxd     ;
extern bool sfc_tau_fxd     ;
extern bool dodamping       ;
extern bool docloud         ;
extern bool docam_sfc_fluxes;
extern bool doprecip        ;
extern bool dosgs           ;
extern bool docoriolis      ;
extern bool dosurface       ;
extern bool dowallx         ;
extern bool dowally         ;
extern bool docolumn        ;
extern bool dotracers       ;
extern bool dosmoke         ;

extern bool advect_sgs;
extern bool dosmagor  ;

extern real vrain;
extern real vsnow;
extern real vgrau;
extern real crain;
extern real csnow;
extern real cgrau;

extern real gam3 ; 
extern real gams1;
extern real gams2;
extern real gams3;
extern real gamg1;
extern real gamg2;
extern real gamg3;
extern real gamr1;
extern real gamr2;
extern real gamr3;

extern real a_bg;
extern real a_pr;
extern real a_gr;

extern bool crm_accel_uv;
extern bool use_crm_accel;
extern real crm_accel_factor;

extern real4d tabs            ;
extern real4d qv              ;
extern real4d qcl             ;
extern real4d qpl             ;
extern real4d qci             ;
extern real4d qpi             ;
extern real5d dudt            ;
extern real5d dvdt            ;
extern real5d dwdt            ;
extern real4d misc            ;
extern real3d fluxbu          ;
extern real3d fluxbv          ;
extern real3d fluxbt          ;
extern real3d fluxbq          ;
extern real3d fluxtu          ;
extern real3d fluxtv          ;
extern real3d fluxtt          ;
extern real3d fluxtq          ;
extern real3d fzero           ;
extern real3d precsfc         ;
extern real3d precssfc        ;
extern real2d t0              ;
extern real2d q0              ;
extern real2d qv0             ;
extern real2d tabs0           ;
extern real2d tv0             ;
extern real2d u0              ;
extern real2d v0              ;
extern real2d tg0             ;
extern real2d qg0             ;
extern real2d ug0             ;
extern real2d vg0             ;
extern real2d p0              ;
extern real2d tke0            ;
extern real2d t01             ;
extern real2d q01             ;
extern real2d qp0             ;
extern real2d qn0             ;
extern real2d prespot         ;
extern real2d rho             ;
extern real2d rhow            ;
extern real2d bet             ;
extern real2d gamaz           ;
extern real2d wsub            ;
extern real2d qtend           ;
extern real2d ttend           ;
extern real2d utend           ;
extern real2d vtend           ;
extern real2d fcorzy          ;
extern real3d latitude        ;
extern real3d longitude       ;
extern real3d prec_xy         ;
extern real3d pw_xy           ;
extern real3d cw_xy           ;
extern real3d iw_xy           ;
extern real3d cld_xy          ;
extern real3d u200_xy         ;
extern real3d usfc_xy         ;
extern real3d v200_xy         ;
extern real3d vsfc_xy         ;
extern real3d w500_xy         ;
extern real2d w_max           ;
extern real2d u_max           ;
extern real2d twsb            ;
extern real2d precflux        ;
extern real2d uwle            ;
extern real2d uwsb            ;
extern real2d vwle            ;
extern real2d vwsb            ;
extern real2d tkelediss       ;
extern real2d tdiff           ;
extern real2d tlat            ;
extern real2d tlatqi          ;
extern real2d qifall          ;
extern real2d qpfall          ;
extern real3d total_water_evap;
extern real3d total_water_prec;
extern real4d CF3D            ;
extern real3d u850_xy         ;
extern real3d v850_xy         ;
extern real3d psfc_xy         ;
extern real3d swvp_xy         ;
extern real3d cloudtopheight  ;
extern real3d echotopheight   ;
extern real3d cloudtoptemp    ;

extern real1d fcorz           ;
extern real1d fcor            ;
extern real1d longitude0      ;
extern real1d latitude0       ;
extern real1d z0              ;
extern real1d uhl             ;
extern real1d vhl             ;
extern real1d taux0           ;
extern real1d tauy0           ;

extern real2d z               ;
extern real2d pres            ;
extern real2d zi              ;
extern real2d presi           ;
extern real2d adz             ;
extern real2d adzw            ;
extern real1d dt3             ;
extern real1d dz              ;

extern real2d grdf_x          ;
extern real2d grdf_y          ;
extern real2d grdf_z          ;
extern real2d tkesbbuoy       ;
extern real2d tkesbshear      ;
extern real2d tkesbdiss       ;

extern real4d fluxbmk         ;
extern real4d fluxtmk         ;
extern real3d mkwle           ;
extern real3d mkwsb           ;
extern real3d mkadv           ;
extern real3d mkdiff          ;
extern real4d qn              ;
extern real2d qpsrc           ;
extern real2d qpevp           ;
extern intHost1d flag_precip  ;
extern int3d flag_top         ;

extern real2d accrsc          ;
extern real2d accrsi          ;
extern real2d accrrc          ;
extern real2d coefice         ;
extern real2d accrgc          ;
extern real2d accrgi          ;
extern real2d evaps1          ;
extern real2d evaps2          ;
extern real2d evapr1          ;
extern real2d evapr2          ;
extern real2d evapg1          ;
extern real2d evapg2          ;

extern real2d t00             ;
extern real2d tln             ;
extern real2d qln             ;
extern real2d qccln           ;
extern real2d qiiln           ;
extern real2d uln             ;
extern real2d vln             ;
extern real3d cwp             ;
extern real3d cwph            ;
extern real3d cwpm            ;
extern real3d cwpl            ;
extern real3d cltemp          ;
extern real3d cmtemp          ;
extern real3d chtemp          ;
extern real3d cttemp          ;
extern real2d dd_crm          ;
extern real2d mui_crm         ;
extern real2d mdi_crm         ;
extern real1d ustar           ;
extern real1d wnd             ;
extern real2d qtot            ;
extern real1d colprec         ;
extern real1d colprecs        ;
extern real1d bflx            ;

extern real1d crm_input_bflxls; 
extern real1d crm_input_wndls ;
extern real2d crm_input_zmid  ;
extern real2d crm_input_zint  ;
extern real2d crm_input_pmid  ;
extern real2d crm_input_pint  ;
extern real2d crm_input_pdel  ;
extern real2d crm_input_ul    ;
extern real2d crm_input_vl    ;
extern real2d crm_input_tl    ;
extern real2d crm_input_qccl  ;
extern real2d crm_input_qiil  ;
extern real2d crm_input_ql    ;
extern real1d crm_input_tau00 ;
extern real4d crm_state_u_wind;
extern real4d crm_state_v_wind;
extern real4d crm_state_w_wind; 
extern real4d crm_state_temperature;
extern real4d crm_state_qt;
extern real4d crm_state_qp;
extern real4d crm_state_qn;
extern real4d crm_rad_qrad;
extern real4d crm_rad_temperature;
extern real4d crm_rad_qv; 
extern real4d crm_rad_qc; 
extern real4d crm_rad_qi; 
extern real4d crm_rad_cld; 
extern real1d crm_output_subcycle_factor;
extern real1d crm_output_prectend;
extern real1d crm_output_precstend; 
extern real2d crm_output_cld; 
extern real2d crm_output_cldtop; 
extern real2d crm_output_gicewp;
extern real2d crm_output_gliqwp; 
extern real2d crm_output_mctot; 
extern real2d crm_output_mcup; 
extern real2d crm_output_mcdn; 
extern real2d crm_output_mcuup; 
extern real2d crm_output_mcudn;
extern real2d crm_output_qc_mean;
extern real2d crm_output_qi_mean;
extern real2d crm_output_qs_mean;
extern real2d crm_output_qg_mean;
extern real2d crm_output_qr_mean;
extern real2d crm_output_mu_crm;
extern real2d crm_output_md_crm;
extern real2d crm_output_eu_crm;
extern real2d crm_output_du_crm;
extern real2d crm_output_ed_crm;
extern real2d crm_output_flux_qt; 
extern real2d crm_output_flux_u;
extern real2d crm_output_flux_v;
extern real2d crm_output_fluxsgs_qt;
extern real2d crm_output_tkez; 
extern real2d crm_output_tkew; 
extern real2d crm_output_tkesgsz; 
extern real2d crm_output_tkz; 
extern real2d crm_output_flux_qp; 
extern real2d crm_output_precflux; 
extern real2d crm_output_qt_trans; 
extern real2d crm_output_qp_trans; 
extern real2d crm_output_qp_fall; 
extern real2d crm_output_qp_evp; 
extern real2d crm_output_qp_src; 
extern real2d crm_output_qt_ls; 
extern real2d crm_output_t_ls; 
extern real1d crm_output_jt_crm; 
extern real1d crm_output_mx_crm; 
extern real1d crm_output_cltot; 
extern real1d crm_output_clhgh; 
extern real1d crm_output_clmed; 
extern real1d crm_output_cllow; 
extern real2d crm_output_sltend; 
extern real2d crm_output_qltend; 
extern real2d crm_output_qcltend; 
extern real2d crm_output_qiltend;
#ifdef MMF_MOMENTUM_FEEDBACK
extern real2d crm_output_ultend; 
extern real2d crm_output_vltend; 
#endif
extern real4d crm_output_tk;
extern real4d crm_output_tkh; 
extern real4d crm_output_qcl; 
extern real4d crm_output_qci; 
extern real4d crm_output_qpl; 
extern real4d crm_output_qpi; 
extern real1d crm_output_z0m; 
extern real1d crm_output_taux; 
extern real1d crm_output_tauy;
extern real1d crm_output_precc; 
extern real1d crm_output_precl;
extern real1d crm_output_precsc; 
extern real1d crm_output_precsl; 
extern real3d crm_output_prec_crm; 
extern real2d crm_clear_rh;
extern int2d  crm_clear_rh_cnt;
extern real1d lat0; 
extern real1d long0;
extern int1d  gcolp;

extern real factor_xy;
extern real factor_xyt;
extern real idt_gl;
extern real crm_nx_rad_fac;
extern real crm_ny_rad_fac;
extern int  ptop;
extern real crm_run_time;
extern real icrm_run_time;
extern real dt_glob;


extern bool crm_accel_ceaseflag;

extern int igstep;


