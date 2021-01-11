#include "pre_timeloop.h"

void pre_timeloop() {
  auto &latitude0                = :: latitude0;
  auto &longitude0               = :: longitude0;
  auto &lat0                     = :: lat0;
  auto &long0                    = :: long0;
  auto &idt_gl                   = :: idt_gl;
  auto &factor_xy                = :: factor_xy;
  auto &crm_rad_temperature      = :: crm_rad_temperature;
  auto &crm_rad_qv               = :: crm_rad_qv;
  auto &crm_rad_qc               = :: crm_rad_qc;
  auto &crm_rad_qi               = :: crm_rad_qi;
  auto &crm_rad_cld              = :: crm_rad_cld;
  auto &crm_clear_rh             = :: crm_clear_rh;
  auto &crm_clear_rh_cnt         = :: crm_clear_rh_cnt;
  auto &bflx                     = :: bflx;
  auto &wnd                      = :: wnd;
  auto &crm_input_bflxls         = :: crm_input_bflxls;
  auto &crm_input_wndls          = :: crm_input_wndls;
  auto &fcor                     = :: fcor;
  auto &fcorz                    = :: fcorz;
  auto &zi                       = :: zi;
  auto &crm_input_zint           = :: crm_input_zint;
  auto &presi                    = :: presi;
  auto &crm_input_pint           = :: crm_input_pint;
  auto &adzw                     = :: adzw;
  auto &fcory                    = :: fcory;
  auto &fcorzy                   = :: fcorzy;
  auto &latitude                 = :: latitude;
  auto &longitude                = :: longitude;
  auto &z                        = :: z;
  auto &crm_input_zmid           = :: crm_input_zmid;
  auto &crm_input_pmid           = :: crm_input_pmid;
  auto &pres                     = :: pres;
  auto &prespot                  = :: prespot;
  auto &bet                      = :: bet;
  auto &crm_input_tl             = :: crm_input_tl;
  auto &gamaz                    = :: gamaz;
  auto &dz                       = :: dz;
  auto &adz                      = :: adz;
  auto &rho                      = :: rho;
  auto &crm_input_pdel           = :: crm_input_pdel;
  auto &u                        = :: u;
  auto &v                        = :: v;
  auto &w                        = :: w;
  auto &tabs                     = :: tabs;
  auto &crm_state_u_wind         = :: crm_state_u_wind;
  auto &crm_state_v_wind         = :: crm_state_v_wind;
  auto &crm_state_w_wind         = :: crm_state_w_wind;
  auto &crm_state_temperature    = :: crm_state_temperature; 
  auto &micro_field              = :: micro_field;
  auto &crm_state_qt             = :: crm_state_qt;
  auto &crm_state_qp             = :: crm_state_qp;
  auto &crm_state_qn             = :: crm_state_qn;
  auto &qn                       = :: qn;
  auto &colprec                  = :: colprec;
  auto &colprecs                 = :: colprecs;
  auto &u0                       = :: u0;
  auto &v0                       = :: v0;
  auto &t0                       = :: t0;
  auto &t00                      = :: t00;
  auto &tabs0                    = :: tabs0;
  auto &q0                       = :: q0;
  auto &qv0                      = :: qv0;
  auto &qn0                      = :: qn0;
  auto &qp0                      = :: qp0;
  auto &tke0                     = :: tke0;
  auto &t                        = :: t;
  auto &qcl                      = :: qcl;
  auto &qci                      = :: qci;
  auto &qpl                      = :: qpl;
  auto &qpi                      = :: qpi;
  auto &sgs_field                = :: sgs_field;
  auto &uln                      = :: uln;
  auto &vln                      = :: vln;
  auto &ttend                    = :: ttend;
  auto &qtend                    = :: qtend;
  auto &crm_input_qccl           = :: crm_input_qccl;
  auto &crm_input_qiil           = :: crm_input_qiil;
  auto &utend                    = :: utend;
  auto &vtend                    = :: vtend;
  auto &ug0                      = :: ug0;
  auto &vg0                      = :: vg0;
  auto &tg0                      = :: tg0;
  auto &qg0                      = :: qg0;
  auto &crm_input_ql             = :: crm_input_ql;
  auto &uhl                      = :: uhl;
  auto &vhl                      = :: vhl;
  auto &crm_input_tau00          = :: crm_input_tau00;
  auto &ustar                    = :: ustar;
  auto &z0                       = :: z0;
  auto &crm_output_subcycle_factor = :: crm_output_subcycle_factor;
  auto &rhow                     = :: rhow;
  auto &qv                       = :: qv;
  auto &crm_output_prectend      = :: crm_output_prectend;
  auto &crm_output_precstend     = :: crm_output_precstend;
  auto &crm_input_ul             = :: crm_input_ul;
  auto &crm_input_vl             = :: crm_input_vl;
  auto &crm_output_cld           = :: crm_output_cld; 
  auto &crm_output_cldtop        = :: crm_output_cldtop; 
  auto &crm_output_gicewp        = :: crm_output_gicewp; 
  auto &crm_output_gliqwp        = :: crm_output_gliqwp; 
  auto &crm_output_mctot         = :: crm_output_mctot; 
  auto &crm_output_mcup          = :: crm_output_mcup; 
  auto &crm_output_mcdn          = :: crm_output_mcdn; 
  auto &crm_output_mcuup         = :: crm_output_mcuup; 
  auto &crm_output_mcudn         = :: crm_output_mcudn; 
  auto &crm_output_qc_mean       = :: crm_output_qc_mean; 
  auto &crm_output_qi_mean       = :: crm_output_qi_mean; 
  auto &crm_output_qs_mean       = :: crm_output_qs_mean; 
  auto &crm_output_qg_mean       = :: crm_output_qg_mean; 
  auto &crm_output_qr_mean       = :: crm_output_qr_mean; 
  auto &crm_output_mu_crm        = :: crm_output_mu_crm; 
  auto &crm_output_md_crm        = :: crm_output_md_crm; 
  auto &crm_output_eu_crm        = :: crm_output_eu_crm; 
  auto &crm_output_du_crm        = :: crm_output_du_crm; 
  auto &crm_output_ed_crm        = :: crm_output_ed_crm; 
  auto &crm_output_flux_qt       = :: crm_output_flux_qt; 
  auto &crm_output_flux_u        = :: crm_output_flux_u; 
  auto &crm_output_flux_v        = :: crm_output_flux_v; 
  auto &crm_output_fluxsgs_qt    = :: crm_output_fluxsgs_qt;
  auto &crm_output_tkez          = :: crm_output_tkez; 
  auto &crm_output_tkew          = :: crm_output_tkew; 
  auto &crm_output_tkesgsz       = :: crm_output_tkesgsz; 
  auto &crm_output_tkz           = :: crm_output_tkz; 
  auto &crm_output_flux_qp       = :: crm_output_flux_qp; 
  auto &crm_output_precflux      = :: crm_output_precflux; 
  auto &crm_output_qt_trans      = :: crm_output_qt_trans; 
  auto &crm_output_qp_trans      = :: crm_output_qp_trans; 
  auto &crm_output_qp_fall       = :: crm_output_qp_fall; 
  auto &crm_output_qp_evp        = :: crm_output_qp_evp; 
  auto &crm_output_qp_src        = :: crm_output_qp_src; 
  auto &crm_output_qt_ls         = :: crm_output_qt_ls; 
  auto &crm_output_t_ls          = :: crm_output_t_ls; 
  auto &dd_crm                   = :: dd_crm; 
  auto &mui_crm                  = :: mui_crm; 
  auto &mdi_crm                  = :: mdi_crm; 
  auto &crm_output_jt_crm        = :: crm_output_jt_crm;
  auto &crm_output_mx_crm        = :: crm_output_mx_crm;
  auto &ncrms                    = :: ncrms;

  crm_accel_ceaseflag = false;

  //Loop over "vector columns"
  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    // latitude0 (icrm) = get_rlat_p(lchnk, icol(icrm)) * 57.296_r8
    // longitude0(icrm) = get_rlon_p(lchnk, icol(icrm)) * 57.296_r8
    latitude0 (icrm) = lat0 (icrm);
    longitude0(icrm) = long0(icrm);
  });

//-----------------------------------------------

  dostatis  = false;    // no statistics are collected.
  idt_gl    = 1.0/dt_glob;
  ptop      = plev-nzm+1;
  factor_xy = 1.0/( (real) nx * (real) ny );

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<crm_ny_rad; j++) {
  //     for (int i=0; i<crm_nx_rad; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,crm_ny_rad,crm_nx_rad,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    crm_rad_temperature(k,j,i,icrm) = 0.0;
    crm_rad_qv(k,j,i,icrm)  = 0.0;
    crm_rad_qc(k,j,i,icrm)  = 0.0;
    crm_rad_qi(k,j,i,icrm)  = 0.0;
    crm_rad_cld(k,j,i,icrm) = 0.0;
  });
  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    bflx(icrm) = crm_input_bflxls(icrm);
    wnd (icrm) = crm_input_wndls (icrm);
  });

//-----------------------------------------

  task_init();
  setparm();

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    fcor(icrm)= 4.0*pi/86400.0*sin(latitude0(icrm)*pi/180.0);
    fcorz(icrm) = sqrt(4.0*pow((2*pi/(3600.*24.)),2)-pow(fcor(icrm),2));
    zi(nz-1,icrm) = crm_input_zint(plev-nz+1,icrm)-crm_input_zint(plev,icrm); //+++mhwang, 2012-02-04
    presi(nz-1,icrm) = crm_input_pint(plev-nz+1,icrm)/100.0;
    adzw(0,icrm) = 1.0;
  });
  // for (int j=0; j<ny+1; j++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(ny+1,ncrms) , YAKL_LAMBDA (int j, int icrm) {
    fcory(j,icrm) = fcor(icrm);
  });
  // for (int j=0; j<ny; j++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(ny,ncrms) , YAKL_LAMBDA (int j, int icrm) {
    fcorzy(j,icrm) = fcorz(icrm);
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    latitude(j,i,icrm) = latitude0(icrm);
    longitude(j,i,icrm) = longitude0(icrm);
  });

  // Create CRM vertical grid and initialize some vertical reference arrays:
  // for (int k=0; k<nzm; k++) {
  //  for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    z(k,icrm) = crm_input_zmid(plev-(k+1),icrm) - crm_input_zint(plev,icrm);
    zi(k,icrm) = crm_input_zint(plev-(k+1)+1,icrm)- crm_input_zint(plev,icrm);
    pres(k,icrm) = crm_input_pmid(plev-(k+1),icrm)/100.0;
    presi(k,icrm) = crm_input_pint(plev-(k+1)+1,icrm)/100.0;
    prespot(k,icrm)=pow((1000.0/pres(k,icrm)),(rgas/cp));
    bet(k,icrm) = ggr/crm_input_tl(plev-(k+1),icrm);
    gamaz(k,icrm)=ggr/cp*z(k,icrm);
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    dz(icrm) = 0.5*(z(0,icrm)+z(1,icrm));
  });

  // for (int k=0; k<nzm-1; k++) {
  //  for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm-1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    k=k+1;
    adzw(k,icrm) = (z(k,icrm)-z(k-1,icrm))/dz(icrm);
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    adzw(nz-1,icrm) = adzw(nzm-1,icrm);
  });
  
  // for (int k=0; k<nzm-1; k++) {
  //  for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    adz(k,icrm)=(zi(k+1,icrm)-zi(k,icrm))/dz(icrm);
    rho(k,icrm) = crm_input_pdel(plev-(k+1),icrm)/ggr/(adz(k,icrm)*dz(icrm));
  });

  // for (int k=0; k<nzm-1; k++) {
  //  for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm-1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    int kInd = k+1;
    rhow(kInd,icrm) = (crm_input_pmid(plev-(k+1),icrm)-crm_input_pmid(plev-(k+1)-1,icrm))/ggr/(adzw(kInd,icrm)*dz(icrm));
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    rhow(0,icrm) = 2.0*rhow(1,icrm) - rhow(2,icrm);
    rhow(nz-1,icrm)= 2.0*rhow(nzm-1,icrm) - rhow(nzm-2,icrm);
  });

  // Initialize clear air relative humidity for aerosol water uptake
  // for (int icrm=0; icrm<ncrms; icrm++) {
  //   for (int k=0; k<nzm; k++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    crm_clear_rh(k,icrm)     = 0.0 ;
    crm_clear_rh_cnt(k,icrm) = 0 ;
  });

  //  Initialize CRM fields:
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    u(k,j+offy_u,i+offx_u,icrm) = crm_state_u_wind(k,j,i,icrm);
    v(k,j+offy_v,i+offx_v,icrm) = crm_state_v_wind(k,j,i,icrm)*YES3D;
    w(k,j+offy_w,i+offx_w,icrm) = crm_state_w_wind(k,j,i,icrm);
    tabs(k,j,i,icrm) = crm_state_temperature(k,j,i,icrm);
  });

  // limit the velocity at the very first step:
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (u(0,offy_u,offx_u,0) == u(0,offy_u,1+offx_u,0) && u(1,offy_u,2+offx_u,0) == u(1,offy_u,3+offx_u,0)) {
      u(k,j+offy_u,i+offx_u,icrm) = min( UMAX, max(-UMAX,u(k,j+offy_u,i+offx_u,icrm)) );
      v(k,j+offy_v,i+offx_v,icrm) = min( UMAX, max(-UMAX,v(k,j+offy_v,i+offx_v,icrm)) )*YES3D;
    }
  });

  // Populate microphysics array from crm_state
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    micro_field(0,k,j+offy_s,i+offx_s,icrm) = crm_state_qt(k,j,i,icrm);
    micro_field(1,k,j+offy_s,i+offx_s,icrm) = crm_state_qp(k,j,i,icrm);
    qn(k,j,i,icrm) = crm_state_qn(k,j,i,icrm);
  });

  micro_init();
  sgs_init();
  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    colprec (icrm)=0.0;
    colprecs(icrm)=0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    u0   (k,icrm)=0.0;
    v0   (k,icrm)=0.0;
    t0   (k,icrm)=0.0;
    t00  (k,icrm)=0.0;
    tabs0(k,icrm)=0.0;
    q0   (k,icrm)=0.0;
    qv0  (k,icrm)=0.0;
    qn0  (k,icrm)=0.0;
    qp0  (k,icrm)=0.0;
    tke0 (k,icrm)=0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    t(k,j+offy_s,i+offx_s,icrm) = tabs(k,j,i,icrm)+gamaz(k,icrm)-fac_cond*qcl(k,j,i,icrm)-fac_sub*qci(k,j,i,icrm) -
                                                                 fac_cond*qpl(k,j,i,icrm)-fac_sub*qpi(k,j,i,icrm);

    real tmp = (qpl(k,j,i,icrm)+qpi(k,j,i,icrm))*crm_input_pdel(plev-(k+1),icrm);
    yakl::atomicAdd(colprec(icrm) , tmp);

    tmp = qpi(k,j,i,icrm)*crm_input_pdel(plev-(k+1),icrm);
    yakl::atomicAdd(colprecs(icrm) , tmp);
    yakl::atomicAdd(u0(k,icrm) , u(k,j+offy_u,i+offx_u,icrm));
    yakl::atomicAdd(v0(k,icrm) , v(k,j+offy_v,i+offx_v,icrm));
    yakl::atomicAdd(t0(k,icrm) , t(k,j+offy_s,i+offx_s,icrm));

    tmp = t(k,j+offy_s,i+offx_s,icrm)+fac_cond*qpl(k,j,i,icrm)+fac_sub*qpi(k,j,i,icrm);
    yakl::atomicAdd(t00(k,icrm) , tmp);
    yakl::atomicAdd(tabs0(k,icrm) , tabs(k,j,i,icrm));

    tmp = qv(k,j,i,icrm)+qcl(k,j,i,icrm)+qci(k,j,i,icrm);
    yakl::atomicAdd(q0(k,icrm) , tmp);
    yakl::atomicAdd(qv0(k,icrm) , qv(k,j,i,icrm));

    tmp = qcl(k,j,i,icrm) + qci(k,j,i,icrm);
    yakl::atomicAdd(qn0(k,icrm) , tmp);

    tmp = qpl(k,j,i,icrm) + qpi(k,j,i,icrm);
    yakl::atomicAdd(qp0(k,icrm) , tmp);
    yakl::atomicAdd(tke0(k,icrm) , sgs_field(0,k,j+offy_s,i+offx_s,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    u0   (k,icrm) = u0   (k,icrm) * factor_xy;
    v0   (k,icrm) = v0   (k,icrm) * factor_xy;
    t0   (k,icrm) = t0   (k,icrm) * factor_xy;
    t00  (k,icrm) = t00  (k,icrm) * factor_xy;
    tabs0(k,icrm) = tabs0(k,icrm) * factor_xy;
    q0   (k,icrm) = q0   (k,icrm) * factor_xy;
    qv0  (k,icrm) = qv0  (k,icrm) * factor_xy;
    qn0  (k,icrm) = qn0  (k,icrm) * factor_xy;
    qp0  (k,icrm) = qp0  (k,icrm) * factor_xy;
    tke0 (k,icrm) = tke0 (k,icrm) * factor_xy;
    int l = plev-(k+1);
    uln  (l,icrm) = min( UMAX, max(-UMAX,crm_input_ul(l,icrm)) );
    vln  (l,icrm) = min( UMAX, max(-UMAX,crm_input_vl(l,icrm)) )*YES3D;
    ttend(k,icrm) = (crm_input_tl(l,icrm)+gamaz(k,icrm)- fac_cond*(crm_input_qccl(l,icrm)+crm_input_qiil(l,icrm))-
                                                         fac_fus*crm_input_qiil(l,icrm)-t00(k,icrm))*idt_gl;
    qtend(k,icrm) = (crm_input_ql(l,icrm)+crm_input_qccl(l,icrm)+crm_input_qiil(l,icrm)-q0(k,icrm))*idt_gl;
    utend(k,icrm) = (uln(l,icrm)-u0(k,icrm))*idt_gl;
    vtend(k,icrm) = (vln(l,icrm)-v0(k,icrm))*idt_gl;
    ug0  (k,icrm) = uln(l,icrm);
    vg0  (k,icrm) = vln(l,icrm);
    tg0  (k,icrm) = crm_input_tl(l,icrm)+gamaz(k,icrm)-fac_cond*crm_input_qccl(l,icrm)-fac_sub*crm_input_qiil(l,icrm);
    qg0  (k,icrm) = crm_input_ql(l,icrm)+crm_input_qccl(l,icrm)+crm_input_qiil(l,icrm);
  });

  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    uhl(icrm) = u0(0,icrm);
    vhl(icrm) = v0(0,icrm);
    // estimate roughness length assuming logarithmic profile of velocity near the surface:
    ustar(icrm) = sqrt(crm_input_tau00(icrm)/rho(0,icrm));
    //z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),ustar(icrm))
    z0_est(z(0,icrm),bflx(icrm),wnd(icrm),ustar(icrm),z0(icrm));
    z0(icrm) = max(0.00001,min(1.0,z0(icrm)));
    crm_output_subcycle_factor(icrm) = 0.0;
    crm_output_prectend (icrm)=colprec (icrm);
    crm_output_precstend(icrm)=colprecs(icrm);
  });

//---------------------------------------------------
  parallel_for( SimpleBounds<2>(plev+1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    if (k < plev) {
      crm_output_cld       (k,icrm) = 0.0;
      crm_output_cldtop    (k,icrm) = 0.0;
      crm_output_gicewp    (k,icrm) = 0.0;
      crm_output_gliqwp    (k,icrm) = 0.0;
      crm_output_mctot     (k,icrm) = 0.0;
      crm_output_mcup      (k,icrm) = 0.0;
      crm_output_mcdn      (k,icrm) = 0.0;
      crm_output_mcuup     (k,icrm) = 0.0;
      crm_output_mcudn     (k,icrm) = 0.0;
      crm_output_qc_mean   (k,icrm) = 0.0;
      crm_output_qi_mean   (k,icrm) = 0.0;
      crm_output_qs_mean   (k,icrm) = 0.0;
      crm_output_qg_mean   (k,icrm) = 0.0;
      crm_output_qr_mean   (k,icrm) = 0.0;
      crm_output_mu_crm    (k,icrm) = 0.0;
      crm_output_md_crm    (k,icrm) = 0.0;
      crm_output_eu_crm    (k,icrm) = 0.0;
      crm_output_du_crm    (k,icrm) = 0.0;
      crm_output_ed_crm    (k,icrm) = 0.0;
      crm_output_flux_qt   (k,icrm) = 0.0;
      crm_output_flux_u    (k,icrm) = 0.0;
      crm_output_flux_v    (k,icrm) = 0.0;
      crm_output_fluxsgs_qt(k,icrm) = 0.0;
      crm_output_tkez      (k,icrm) = 0.0;
      crm_output_tkew      (k,icrm) = 0.0;
      crm_output_tkesgsz   (k,icrm) = 0.0;
      crm_output_tkz       (k,icrm) = 0.0;
      crm_output_flux_qp   (k,icrm) = 0.0;
      crm_output_precflux  (k,icrm) = 0.0;
      crm_output_qt_trans  (k,icrm) = 0.0;
      crm_output_qp_trans  (k,icrm) = 0.0;
      crm_output_qp_fall   (k,icrm) = 0.0;
      crm_output_qp_evp    (k,icrm) = 0.0;
      crm_output_qp_src    (k,icrm) = 0.0;
      crm_output_qt_ls     (k,icrm) = 0.0;
      crm_output_t_ls      (k,icrm) = 0.0;
      dd_crm               (k,icrm) = 0.0;
    }
    mui_crm(k,icrm) = 0.0;
    mdi_crm(k,icrm) = 0.0;
  });

  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    crm_output_jt_crm(icrm) = 0.0;
    crm_output_mx_crm(icrm) = 0.0;
  });

//--------------------------------------------------
  if (doprecip) { 
    precip_init();
  }
  
  if ( igstep <= 1 ) { setperturb(); }

  if ( nx%crm_nx_rad==0 || ny%crm_ny_rad==0  ) {
    crm_nx_rad_fac = static_cast<real>(crm_nx_rad)/static_cast<real>(nx);
    crm_ny_rad_fac = static_cast<real>(crm_ny_rad)/static_cast<real>(ny);
  } else {
    std::cout << "crm_nx_rad and crm_ny_rad need to be divisible by nx and ny";
    exit(-1);
  }

  //perturb_arrays()

  nstop = dt_glob/dt;
  dt = dt_glob/((real) nstop);

  crm_run_time  = dt_glob;
  icrm_run_time = 1.0/crm_run_time;

  if (use_crm_accel) {
    crm_accel_nstop(nstop);  // reduce nstop by factor of (1 + crm_accel_factor)
  }

}
