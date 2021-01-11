
#include "post_timeloop.h"

void post_timeloop() {
  auto &crm_rad_temperature     = :: crm_rad_temperature;
  auto &crm_rad_qv              = :: crm_rad_qv;
  auto &crm_rad_qc              = :: crm_rad_qc;
  auto &crm_rad_qi              = :: crm_rad_qi;
  auto &crm_rad_cld             = :: crm_rad_cld;
  auto &tln                     = :: tln;
  auto &qln                     = :: qln;
  auto &qccln                   = :: qccln;
  auto &qiiln                   = :: qiiln;
  auto &uln                     = :: uln;
  auto &vln                     = :: vln;
  auto &crm_input_tl            = :: crm_input_tl;
  auto &crm_input_ql            = :: crm_input_ql;
  auto &crm_input_qccl          = :: crm_input_qccl;
  auto &crm_input_qiil          = :: crm_input_qiil;
  auto &crm_input_ul            = :: crm_input_ul;
  auto &crm_input_vl            = :: crm_input_vl;
  auto &colprec                 = :: colprec;
  auto &colprecs                = :: colprecs;
  auto &qpl                     = :: qpl;
  auto &qpi                     = :: qpi;
  auto &crm_input_pdel          = :: crm_input_pdel;
  auto &tabs                    = :: tabs;
  auto &qv                      = :: qv;
  auto &u                       = :: u;
  auto &v                       = :: v;
  auto &qcl                     = :: qcl;
  auto &qci                     = :: qci;
  auto &factor_xyt              = :: factor_xyt;
  auto &factor_xy               = :: factor_xy;
  auto &nstop                   = :: nstop;
  auto &crm_output_sltend       = :: crm_output_sltend;
  auto &crm_output_qltend       = :: crm_output_qltend;
  auto &crm_output_qcltend      = :: crm_output_qcltend;
  auto &crm_output_qiltend      = :: crm_output_qiltend;
#ifdef MMF_MOMENTUM_FEEDBACK
  auto &crm_output_ultend       = :: crm_output_ultend;
  auto &crm_output_vltend       = :: crm_output_vltend;
#endif
  auto &icrm_run_time           = :: icrm_run_time;
  auto &crm_output_prectend     = :: crm_output_prectend;
  auto &crm_output_precstend    = :: crm_output_precstend;
  auto &w                       = :: w;
  auto &crm_state_u_wind        = :: crm_state_u_wind;
  auto &crm_state_v_wind        = :: crm_state_v_wind;
  auto &crm_state_w_wind        = :: crm_state_w_wind;
  auto &crm_state_temperature   = :: crm_state_temperature;
  auto &crm_state_qt            = :: crm_state_qt;
  auto &crm_state_qp            = :: crm_state_qp;
  auto &crm_state_qn            = :: crm_state_qn;
  auto &micro_field             = :: micro_field;
  auto &sgs_field_diag          = :: sgs_field_diag;
  auto &qn                      = :: qn;
  auto &crm_output_qcl          = :: crm_output_qcl;
  auto &crm_output_qci          = :: crm_output_qci;
  auto &crm_output_qpl          = :: crm_output_qpl;
  auto &crm_output_qpi          = :: crm_output_qpi;
  auto &crm_output_tk           = :: crm_output_tk;
  auto &crm_output_tkh          = :: crm_output_tkh;
  auto &crm_output_z0m          = :: crm_output_z0m;
  auto &crm_output_taux         = :: crm_output_taux;
  auto &crm_output_tauy         = :: crm_output_tauy;
  auto &z0                      = :: z0;
  auto &taux0                   = :: taux0;
  auto &tauy0                   = :: tauy0;
  auto &crm_output_qc_mean      = :: crm_output_qc_mean;
  auto &crm_output_qi_mean      = :: crm_output_qi_mean;
  auto &crm_output_qr_mean      = :: crm_output_qr_mean;
  auto &crm_output_qg_mean      = :: crm_output_qg_mean;
  auto &crm_output_qs_mean      = :: crm_output_qs_mean;
  auto &crm_output_cld          = :: crm_output_cld;
  auto &crm_output_cldtop       = :: crm_output_cldtop;
  auto &crm_output_gicewp       = :: crm_output_gicewp;
  auto &crm_output_gliqwp       = :: crm_output_gliqwp;
  auto &crm_output_mcup         = :: crm_output_mcup;
  auto &crm_output_mcdn         = :: crm_output_mcdn;
  auto &crm_output_mcuup        = :: crm_output_mcuup;
  auto &crm_output_mcudn        = :: crm_output_mcudn;
  auto &crm_output_mctot        = :: crm_output_mctot;
  auto &crm_output_precc        = :: crm_output_precc;
  auto &crm_output_precl        = :: crm_output_precl;
  auto &crm_output_precsc       = :: crm_output_precsc;
  auto &crm_output_precsl       = :: crm_output_precsl;
  auto &precsfc                 = :: precsfc;
  auto &precssfc                = :: precssfc;
  auto &crm_output_prec_crm     = :: crm_output_prec_crm;
  auto &crm_clear_rh            = :: crm_clear_rh;
  auto &crm_clear_rh_cnt        = :: crm_clear_rh_cnt;
  auto &crm_output_cltot        = :: crm_output_cltot;
  auto &crm_output_clhgh        = :: crm_output_clhgh;
  auto &crm_output_clmed        = :: crm_output_clmed;
  auto &crm_output_cllow        = :: crm_output_cllow;
  auto &crm_output_jt_crm       = :: crm_output_jt_crm;
  auto &crm_output_mx_crm       = :: crm_output_mx_crm;
  auto &crm_output_mu_crm       = :: crm_output_mu_crm;
  auto &crm_output_md_crm       = :: crm_output_md_crm;
  auto &crm_output_eu_crm       = :: crm_output_eu_crm;
  auto &crm_output_du_crm       = :: crm_output_du_crm;
  auto &crm_output_ed_crm       = :: crm_output_ed_crm;
  auto &mui_crm                 = :: mui_crm;
  auto &mdi_crm                 = :: mdi_crm;
  auto &dd_crm                  = :: dd_crm;
  auto &u0                      = :: u0;
  auto &v0                      = :: v0;
  auto &mkwsb                   = :: mkwsb;
  auto &mkwle                   = :: mkwle;
  auto &mkadv                   = :: mkadv;
  auto &mkdiff                  = :: mkdiff;
  auto &qpsrc                   = :: qpsrc;
  auto &qpevp                   = :: qpevp;
  auto &qpfall                  = :: qpfall;
  auto &precflux                = :: precflux;
  auto &crm_output_flux_u       = :: crm_output_flux_u;
  auto &crm_output_flux_v       = :: crm_output_flux_v;
  auto &crm_output_flux_qt      = :: crm_output_flux_qt;
  auto &crm_output_fluxsgs_qt   = :: crm_output_fluxsgs_qt;
  auto &crm_output_flux_qp      = :: crm_output_flux_qp;
  auto &crm_output_qt_trans     = :: crm_output_qt_trans;
  auto &crm_output_qp_trans     = :: crm_output_qp_trans;
  auto &uwle                    = :: uwle;
  auto &vwle                    = :: vwle;
  auto &uwsb                    = :: uwsb;
  auto &vwsb                    = :: vwsb;
  auto &crm_output_tkesgsz      = :: crm_output_tkesgsz;
  auto &crm_output_tkez         = :: crm_output_tkez;
  auto &crm_output_tkew         = :: crm_output_tkew;
  auto &crm_output_tkz          = :: crm_output_tkz;
  auto &crm_output_precflux     = :: crm_output_precflux;
  auto &crm_output_qp_fall      = :: crm_output_qp_fall;
  auto &crm_output_qp_evp       = :: crm_output_qp_evp;
  auto &crm_output_qp_src       = :: crm_output_qp_src;
  auto &crm_output_qt_ls        = :: crm_output_qt_ls;
  auto &crm_output_t_ls         = :: crm_output_t_ls;
  auto &qtend                   = :: qtend;
  auto &ttend                   = :: ttend;
  auto &a_gr                    = :: a_gr;
  auto &dz                      = :: dz;
  auto &dt                      = :: dt;
  auto &ptop                    = :: ptop;
  auto &rhow                    = :: rhow;
  auto &rho                     = :: rho;
  auto &sgs_field               = :: sgs_field;
  auto &dtn                     = :: dtn;
  auto &crm_output_subcycle_factor= :: crm_output_subcycle_factor;
  auto &ncrms                   = :: ncrms;

  factor_xyt = factor_xy/((real) nstop);
  real tmp1 = crm_nx_rad_fac*crm_ny_rad_fac/((real) nstop);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<crm_ny_rad; j++) {
  //     for (int i=0; i<crm_nx_rad; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,crm_ny_rad,crm_nx_rad,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    crm_rad_temperature(k,j,i,icrm) = crm_rad_temperature(k,j,i,icrm) * tmp1;
    crm_rad_qv         (k,j,i,icrm) = crm_rad_qv         (k,j,i,icrm) * tmp1;
    crm_rad_qc         (k,j,i,icrm) = crm_rad_qc         (k,j,i,icrm) * tmp1;
    crm_rad_qi         (k,j,i,icrm) = crm_rad_qi         (k,j,i,icrm) * tmp1;
    crm_rad_cld        (k,j,i,icrm) = crm_rad_cld        (k,j,i,icrm) * tmp1;
  });

  // Convert clear RH sum to average
  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    if (crm_clear_rh_cnt(k,icrm)>0) {
      crm_clear_rh(k,icrm) = crm_clear_rh(k,icrm) / crm_clear_rh_cnt(k,icrm);
    }
  });

  // no CRM tendencies above its top
  // for (int k=0; k<ptop-1; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(ptop-1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    tln  (k,icrm) = crm_input_tl  (k,icrm);
    qln  (k,icrm) = crm_input_ql  (k,icrm);
    qccln(k,icrm) = crm_input_qccl(k,icrm);
    qiiln(k,icrm) = crm_input_qiil(k,icrm);
    uln  (k,icrm) = crm_input_ul  (k,icrm);
    vln  (k,icrm) = crm_input_vl  (k,icrm);
  });

  //  Compute tendencies due to CRM:
  // for (int k=0; k<plev-ptop+1; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(plev-ptop+1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    k = ptop+k-1;
    tln  (k,icrm) = 0.0;
    qln  (k,icrm) = 0.0;
    qccln(k,icrm) = 0.0;
    qiiln(k,icrm) = 0.0;
    uln  (k,icrm) = 0.0;
    vln  (k,icrm) = 0.0;
  });
  
  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    colprec (icrm)=0;
    colprecs(icrm)=0;
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int j=0; j<ny; j++) {
  //      for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int l = plev-(k+1);

    real tmp = (qpl(k,j,i,icrm)+qpi(k,j,i,icrm))*crm_input_pdel(l,icrm);
    yakl::atomicAdd(colprec (icrm) , tmp);

    tmp = qpi(k,j,i,icrm)*crm_input_pdel(l,icrm);
    yakl::atomicAdd(colprecs(icrm) , tmp);
    yakl::atomicAdd(tln(l,icrm) , tabs(k,j,i,icrm));
    yakl::atomicAdd(qln(l,icrm) , qv(k,j,i,icrm));
    yakl::atomicAdd(qccln(l,icrm) , qcl(k,j,i,icrm));
    yakl::atomicAdd(qiiln(l,icrm) , qci(k,j,i,icrm));
    yakl::atomicAdd(uln(l,icrm) , u(k,j+offy_u,i+offx_u,icrm));
    yakl::atomicAdd(vln(l,icrm) , v(k,j+offy_v,i+offx_v,icrm));
  });

  // for (int k=0; k<plev-ptop+1; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(plev-ptop+1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    k = ptop+k-1;

    tln  (k,icrm) = tln  (k,icrm) * factor_xy;
    qln  (k,icrm) = qln  (k,icrm) * factor_xy;
    qccln(k,icrm) = qccln(k,icrm) * factor_xy;
    qiiln(k,icrm) = qiiln(k,icrm) * factor_xy;
    uln  (k,icrm) = uln  (k,icrm) * factor_xy;
    vln  (k,icrm) = vln  (k,icrm) * factor_xy;
  });

  // for (int k=0; k<plev; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(plev,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    crm_output_sltend (k,icrm) = cp * (tln  (k,icrm) - crm_input_tl  (k,icrm)) * icrm_run_time;
    crm_output_qltend (k,icrm) =      (qln  (k,icrm) - crm_input_ql  (k,icrm)) * icrm_run_time;
    crm_output_qcltend(k,icrm) =      (qccln(k,icrm) - crm_input_qccl(k,icrm)) * icrm_run_time;
    crm_output_qiltend(k,icrm) =      (qiiln(k,icrm) - crm_input_qiil(k,icrm)) * icrm_run_time;
#ifdef MMF_MOMENTUM_FEEDBACK
    crm_output_ultend (k,icrm) =      (uln  (k,icrm) - crm_input_ul  (k,icrm)) * icrm_run_time;
    crm_output_vltend (k,icrm) =      (vln  (k,icrm) - crm_input_vl  (k,icrm)) * icrm_run_time;
#endif
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    crm_output_prectend (icrm) = (colprec (icrm)-crm_output_prectend (icrm))/ggr*factor_xy * icrm_run_time;
    crm_output_precstend(icrm) = (colprecs(icrm)-crm_output_precstend(icrm))/ggr*factor_xy * icrm_run_time;
  });

  // don't use CRM tendencies from two crm top levels
  // radiation tendencies are added back after the CRM call (see crm_physics_tend)

  // for (int k=0; k<2; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(2,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    k = ptop+k-1;

    crm_output_sltend (k,icrm) = 0.0;
    crm_output_qltend (k,icrm) = 0.0;
    crm_output_qcltend(k,icrm) = 0.0;
    crm_output_qiltend(k,icrm) = 0.0;
#ifdef MMF_MOMENTUM_FEEDBACK
    crm_output_ultend (k,icrm) = 0.0;
    crm_output_vltend (k,icrm) = 0.0;
#endif
  });

  // Save the last step to the permanent core:

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    crm_state_u_wind(k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm);
    crm_state_v_wind(k,j,i,icrm) = v(k,j+offy_v,i+offx_v,icrm);
    crm_state_w_wind(k,j,i,icrm) = w(k,j+offy_w,i+offx_w,icrm);
    crm_state_temperature(k,j,i,icrm) = tabs(k,j,i,icrm);
    crm_state_qt(k,j,i,icrm) = micro_field(0,k,j+offy_s,i+offx_s,icrm);
    crm_state_qp(k,j,i,icrm) = micro_field(1,k,j+offy_s,i+offx_s,icrm);
    crm_state_qn(k,j,i,icrm) = qn(k,j,i,icrm);
    crm_output_tk(k,j,i,icrm) = sgs_field_diag(0,k,j+offy_d,i+offx_d,icrm);
    crm_output_tkh(k,j,i,icrm) = sgs_field_diag(1,k,j+offy_d,i+offx_d,icrm);
    crm_output_qcl(k,j,i,icrm) = qcl(k,j,i,icrm);
    crm_output_qci(k,j,i,icrm) = qci(k,j,i,icrm);
    crm_output_qpl(k,j,i,icrm) = qpl(k,j,i,icrm);
    crm_output_qpi(k,j,i,icrm) = qpi(k,j,i,icrm);
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    crm_output_z0m (icrm) = z0(icrm);
    crm_output_taux(icrm) = taux0(icrm) / nstop;
    crm_output_tauy(icrm) = tauy0(icrm) / nstop;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int l = plev-(k+1);
    yakl::atomicAdd(crm_output_qc_mean(l,icrm) , qcl(k,j,i,icrm));
    yakl::atomicAdd(crm_output_qi_mean(l,icrm) , qci(k,j,i,icrm));
    yakl::atomicAdd(crm_output_qr_mean(l,icrm) , qpl(k,j,i,icrm));
    real omg = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tgrmin)*a_gr));

    real tmp = qpi(k,j,i,icrm)*omg;
    yakl::atomicAdd(crm_output_qg_mean(l,icrm) , tmp);

    tmp = qpi(k,j,i,icrm)*(1.0-omg);
    yakl::atomicAdd(crm_output_qs_mean(l,icrm) , tmp);
  });

  // for (int k=0; k<plev; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(plev,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    crm_output_cld   (k,icrm) = min( 1.0, crm_output_cld   (k,icrm) * factor_xyt );
    crm_output_cldtop(k,icrm) = min( 1.0, crm_output_cldtop(k,icrm) * factor_xyt );
    crm_output_gicewp(k,icrm) = crm_output_gicewp(k,icrm)*crm_input_pdel(k,icrm)*1000.0/ggr * factor_xyt;
    crm_output_gliqwp(k,icrm) = crm_output_gliqwp(k,icrm)*crm_input_pdel(k,icrm)*1000.0/ggr * factor_xyt;
    crm_output_mcup  (k,icrm) = crm_output_mcup (k,icrm) * factor_xyt;
    crm_output_mcdn  (k,icrm) = crm_output_mcdn (k,icrm) * factor_xyt;
    crm_output_mcuup (k,icrm) = crm_output_mcuup(k,icrm) * factor_xyt;
    crm_output_mcudn (k,icrm) = crm_output_mcudn(k,icrm) * factor_xyt;
    crm_output_mctot (k,icrm) = crm_output_mcup(k,icrm) + crm_output_mcdn(k,icrm) + 
                                crm_output_mcuup(k,icrm) + crm_output_mcudn(k,icrm);

    crm_output_qc_mean(k,icrm) = crm_output_qc_mean(k,icrm) * factor_xy;
    crm_output_qi_mean(k,icrm) = crm_output_qi_mean(k,icrm) * factor_xy;
    crm_output_qs_mean(k,icrm) = crm_output_qs_mean(k,icrm) * factor_xy;
    crm_output_qg_mean(k,icrm) = crm_output_qg_mean(k,icrm) * factor_xy;
    crm_output_qr_mean(k,icrm) = crm_output_qr_mean(k,icrm) * factor_xy;
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    crm_output_precc (icrm) = 0.0;
    crm_output_precl (icrm) = 0.0;
    crm_output_precsc(icrm) = 0.0;
    crm_output_precsl(icrm) = 0.0;
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    precsfc(j,i,icrm) = precsfc(j,i,icrm)*dz(icrm)/dt/((real) nstop);
    precssfc(j,i,icrm) = precssfc(j,i,icrm)*dz(icrm)/dt/((real) nstop);
    if (precsfc(j,i,icrm) > 10.0/86400.0) {
      yakl::atomicAdd(crm_output_precc (icrm) , precsfc (j,i,icrm));
      yakl::atomicAdd(crm_output_precsc(icrm) , precssfc(j,i,icrm));
    } else {
      yakl::atomicAdd(crm_output_precl (icrm) , precsfc (j,i,icrm));
      yakl::atomicAdd(crm_output_precsl(icrm) , precssfc(j,i,icrm));
    }
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    crm_output_prec_crm(j,i,icrm) = precsfc(j,i,icrm)/1000.0;           //mm/s --> m/s
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    crm_output_precc (icrm) = crm_output_precc (icrm)*factor_xy/1000.0;
    crm_output_precl (icrm) = crm_output_precl (icrm)*factor_xy/1000.0;
    crm_output_precsc(icrm) = crm_output_precsc(icrm)*factor_xy/1000.0;
    crm_output_precsl(icrm) = crm_output_precsl(icrm)*factor_xy/1000.0;

    crm_output_cltot(icrm) = crm_output_cltot(icrm) * factor_xyt;
    crm_output_clhgh(icrm) = crm_output_clhgh(icrm) * factor_xyt;
    crm_output_clmed(icrm) = crm_output_clmed(icrm) * factor_xyt;
    crm_output_cllow(icrm) = crm_output_cllow(icrm) * factor_xyt;

    crm_output_jt_crm(icrm) = plev * 1.0;
    crm_output_mx_crm(icrm) = 1.0;
  });

  // for (int k=0; k<plev; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(plev,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    crm_output_mu_crm(k,icrm)=0.5*(mui_crm(k,icrm)+mui_crm(k+1,icrm));
    crm_output_md_crm(k,icrm)=0.5*(mdi_crm(k,icrm)+mdi_crm(k+1,icrm));
    crm_output_mu_crm(k,icrm)=crm_output_mu_crm(k,icrm)*ggr/100.0;          //kg/m2/s --> mb/s
    crm_output_md_crm(k,icrm)=crm_output_md_crm(k,icrm)*ggr/100.0;          //kg/m2/s --> mb/s
    crm_output_eu_crm(k,icrm) = 0.0;
    if (mui_crm(k,icrm)-mui_crm(k+1,icrm) > 0) {
      crm_output_eu_crm(k,icrm)=(mui_crm(k,icrm)-mui_crm(k+1,icrm))*ggr/crm_input_pdel(k,icrm);    // /s
    } else {
      crm_output_du_crm(k,icrm)=-1.0*(mui_crm(k,icrm)-mui_crm(k+1,icrm))*ggr/crm_input_pdel(k,icrm);   // /s
    }
    if (mdi_crm(k+1,icrm)-mdi_crm(k,icrm) < 0) {
      crm_output_ed_crm(k,icrm)=(mdi_crm(k,icrm)-mdi_crm(k+1,icrm))*ggr/crm_input_pdel(k,icrm); // /s
    } else {
      dd_crm(k,icrm)=-1.0*(mdi_crm(k,icrm)-mdi_crm(k+1,icrm))*ggr/crm_input_pdel(k,icrm);   // /s
    }

    real tmp;
    if (abs(crm_output_mu_crm(k,icrm)) > 1.0e-15 || abs(crm_output_md_crm(k,icrm)) > 1.0e-15) {
      tmp = k+1;
      yakl::atomicMin(crm_output_jt_crm(icrm),tmp);

      tmp = k+1;
      yakl::atomicMax(crm_output_mx_crm(icrm),tmp);
    }
  });

  //       Fluxes and other stat:
  
  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    real u2z = 0.0;
    real v2z = 0.0;
    real w2z = 0.0;
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        real tmp1 = (u(k,j+offy_u,i+offx_u,icrm)-u0(k,icrm));
        real tmp2 = (v(k,j+offy_v,i+offx_v,icrm)-v0(k,icrm));
        real tmp3 = w(k+1,j+offy_w,i+offx_w,icrm);
        real tmp4 = w(k,j+offy_w,i+offx_w,icrm);
        u2z = u2z+tmp1*tmp1;
        v2z = v2z+tmp2*tmp2;
        w2z = w2z+0.5*(tmp3*tmp3+tmp4*tmp4);
      }
    }

    real tmp1 = dz(icrm)/rhow(k,icrm);
    real tmp2 = tmp1/dtn; // dtn is calculated inside of the icyc loop. It seems wrong to use it here ???? +++mhwang

    for (int l=0; l<nmicro_fields; l++) {                                           
      mkwsb(l,k,icrm) = mkwsb(l,k,icrm) * tmp1*rhow(k,icrm) * factor_xy/((real) nstop);     //kg/m3/s --> kg/m2/s
      mkwle(l,k,icrm) = mkwle(l,k,icrm) * tmp2*rhow(k,icrm) * factor_xy/((real) nstop);     //kg/m3   --> kg/m2/s
      mkadv(l,k,icrm) = mkadv(l,k,icrm) * factor_xy*icrm_run_time;     // kg/kg  --> kg/kg/s
      mkdiff(l,k,icrm) = mkdiff(l,k,icrm) * factor_xy*icrm_run_time;   // kg/kg  --> kg/kg/s
    }

    // qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux.
    qpsrc(k,icrm) = qpsrc(k,icrm) * factor_xy*icrm_run_time;
    qpevp(k,icrm) = qpevp(k,icrm) * factor_xy*icrm_run_time;
    qpfall(k,icrm) = qpfall(k,icrm) * factor_xy*icrm_run_time;   // kg/kg in M2005 ---> kg/kg/s
    precflux(k,icrm) = precflux(k,icrm) * factor_xy*dz(icrm)/dt/((real) nstop);  //kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/((real) nstop))

    int l = plev-(k+1);
    crm_output_flux_u    (l,icrm) = (uwle(k,icrm) + uwsb(k,icrm))*tmp1*factor_xy/((real) nstop);
    crm_output_flux_v    (l,icrm) = (vwle(k,icrm) + vwsb(k,icrm))*tmp1*factor_xy/((real) nstop);
    crm_output_flux_qt   (l,icrm) = mkwle(0,k,icrm) + mkwsb(0,k,icrm);
    crm_output_fluxsgs_qt(l,icrm) = mkwsb(0,k,icrm);
    crm_output_flux_qp   (l,icrm) = mkwle(1,k,icrm) + mkwsb(1,k,icrm);
    crm_output_qt_trans  (l,icrm) = mkadv(0,k,icrm) + mkdiff(0,k,icrm);
    crm_output_qp_trans  (l,icrm) = mkadv(1,k,icrm) + mkdiff(1,k,icrm);
    real tmp = 0.0;
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        tmp = tmp + sgs_field(0,k,j+offy_s,i+offx_s,icrm);
      }
    }
    crm_output_tkesgsz   (l,icrm)= rho(k,icrm)*tmp*factor_xy;
    crm_output_tkez      (l,icrm)= rho(k,icrm)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + crm_output_tkesgsz(l,icrm);
    crm_output_tkew      (l,icrm)= rho(k,icrm)*0.5*w2z*factor_xy;

    tmp = 0.0;
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        tmp = tmp + sgs_field_diag(0,k,j+offy_d,i+offx_d,icrm);
      }
    }
    crm_output_tkz       (l,icrm) = tmp * factor_xy;
    crm_output_precflux  (l,icrm) = precflux(k,icrm)/1000.0;      // mm/s  -->m/s

    crm_output_qp_fall   (l,icrm) = qpfall(k,icrm);
    crm_output_qp_evp    (l,icrm) = qpevp(k,icrm);
    crm_output_qp_src    (l,icrm) = qpsrc(k,icrm);

    crm_output_qt_ls     (l,icrm) = qtend(k,icrm);
    crm_output_t_ls      (l,icrm) = ttend(k,icrm);
  });

  parallel_for( ncrms , YAKL_LAMBDA(int icrm) {
    crm_output_subcycle_factor(icrm) = crm_output_subcycle_factor(icrm)/((real) nstop);
  });
}


