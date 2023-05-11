#include "post_icycle.h"

void post_icycle() {
  YAKL_SCOPE( cwp                 , :: cwp );
  YAKL_SCOPE( cwph                , :: cwph );
  YAKL_SCOPE( cwpm                , :: cwpm );
  YAKL_SCOPE( cwpl                , :: cwpl );
  YAKL_SCOPE( flag_top            , :: flag_top );
  YAKL_SCOPE( cltemp              , :: cltemp );
  YAKL_SCOPE( cmtemp              , :: cmtemp );
  YAKL_SCOPE( chtemp              , :: chtemp );
  YAKL_SCOPE( cttemp              , :: cttemp );
  YAKL_SCOPE( rho                 , :: rho );
  YAKL_SCOPE( adz                 , :: adz );
  YAKL_SCOPE( dz                  , :: dz );
  YAKL_SCOPE( qcl                 , :: qcl );
  YAKL_SCOPE( qci                 , :: qci );
  YAKL_SCOPE( CF3D                , :: CF3D );
  YAKL_SCOPE( crm_output_cldtop   , :: crm_output_cldtop );
  YAKL_SCOPE( w                   , :: w );
  YAKL_SCOPE( crm_output_mcup     , :: crm_output_mcup );
  YAKL_SCOPE( crm_output_mcuup    , :: crm_output_mcuup );
  YAKL_SCOPE( crm_output_mcdn     , :: crm_output_mcdn );
  YAKL_SCOPE( crm_output_mcudn    , :: crm_output_mcudn );
  YAKL_SCOPE( crm_output_gliqwp   , :: crm_output_gliqwp );
  YAKL_SCOPE( crm_output_gicewp   , :: crm_output_gicewp );
  YAKL_SCOPE( crm_output_cld      , :: crm_output_cld );
  YAKL_SCOPE( crm_rad_temperature , :: crm_rad_temperature );
  YAKL_SCOPE( crm_rad_qv          , :: crm_rad_qv );
  YAKL_SCOPE( crm_rad_qc          , :: crm_rad_qc );
  YAKL_SCOPE( crm_rad_qi          , :: crm_rad_qi );
  YAKL_SCOPE( tabs                , :: tabs );
  YAKL_SCOPE( crm_rad_cld         , :: crm_rad_cld );
  YAKL_SCOPE( crm_clear_rh        , :: crm_clear_rh );
  YAKL_SCOPE( crm_clear_rh_cnt    , :: crm_clear_rh_cnt );
  YAKL_SCOPE( pres                , :: pres );
  YAKL_SCOPE( rhow                , :: rhow );
  YAKL_SCOPE( mui_crm             , :: mui_crm );
  YAKL_SCOPE( mdi_crm             , :: mdi_crm );
  YAKL_SCOPE( crm_output_cltot    , :: crm_output_cltot );
  YAKL_SCOPE( crm_output_clhgh    , :: crm_output_clhgh );
  YAKL_SCOPE( crm_output_clmed    , :: crm_output_clmed );
  YAKL_SCOPE( crm_output_cllow    , :: crm_output_cllow );
  YAKL_SCOPE( qv                  , :: qv );
  YAKL_SCOPE( qpi                 , :: qpi );
  YAKL_SCOPE( qpl                 , :: qpl );
  YAKL_SCOPE( ncrms               , :: ncrms );

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    cwp     (j,i,icrm) = 0.0;
    cwph    (j,i,icrm) = 0.0;
    cwpm    (j,i,icrm) = 0.0;
    cwpl    (j,i,icrm) = 0.0;
    flag_top(j,i,icrm) = 1;
    cltemp  (j,i,icrm) = 0.0;
    cmtemp  (j,i,icrm) = 0.0;
    chtemp  (j,i,icrm) = 0.0;
    cttemp  (j,i,icrm) = 0.0;
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    for (int k=0; k<nzm; k++) {
      int l = plev-(k+1);
      real tmp1 = rho(nz-(k+1)-1,icrm)*adz(nz-(k+1)-1,icrm)*dz(icrm)*(qcl(nz-(k+1)-1,j,i,icrm)+qci(nz-(k+1)-1,j,i,icrm));
      cwp   (j,i,icrm) = cwp(j,i,icrm)+tmp1;
      cttemp(j,i,icrm) = max(CF3D(nz-(k+1)-1,j,i,icrm), cttemp(j,i,icrm));
      if (cwp(j,i,icrm) > cwp_threshold && flag_top(j,i,icrm) == 1) {
        yakl::atomicAdd(crm_output_cldtop(l,icrm), 1.0);
        flag_top(j,i,icrm) = 0;
      }
      if (pres(nz-(k+1)-1,icrm) >= 700.0) {
        cwpl(j,i,icrm) = cwpl(j,i,icrm)+tmp1;
        cltemp(j,i,icrm) = max(CF3D(nz-(k+1)-1,j,i,icrm), cltemp(j,i,icrm));
      } else if (pres(nz-(k+1)-1,icrm) < 400.0) {
        cwph(j,i,icrm) = cwph(j,i,icrm)+tmp1;
        chtemp(j,i,icrm) = max(CF3D(nz-(k+1)-1,j,i,icrm), chtemp(j,i,icrm));
      } else {
        cwpm(j,i,icrm) = cwpm(j,i,icrm)+tmp1;
        cmtemp(j,i,icrm) = max(CF3D(nz-(k+1)-1,j,i,icrm), cmtemp(j,i,icrm));
      }
      tmp1 = rho(k,icrm)*adz(k,icrm)*dz(icrm);
      real tmp;
      if(tmp1*(qcl(k,j,i,icrm)+qci(k,j,i,icrm)) > cwp_threshold) {
         yakl::atomicAdd(crm_output_cld(l,icrm), CF3D(k,j,i,icrm));
         if(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm) > 2*wmin) {
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm)) * CF3D(k,j,i,icrm);
           yakl::atomicAdd(crm_output_mcup(l,icrm), tmp);
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm)) * (1.0 - CF3D(k,j,i,icrm));
           yakl::atomicAdd(crm_output_mcuup(l,icrm) , tmp);
         }
         if(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm) < -2*wmin) {
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm)) * CF3D(k,j,i,icrm);
           yakl::atomicAdd(crm_output_mcdn (l,icrm) , tmp);
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm)) * (1. - CF3D(k,j,i,icrm));
           yakl::atomicAdd(crm_output_mcudn(l,icrm) , tmp);
         }
      } else {
         if(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm) > 2*wmin) {
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm));
           yakl::atomicAdd(crm_output_mcuup(l,icrm) , tmp);
         }
         if(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm) < -2*wmin) {
           tmp = rho(k,icrm)*0.5*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm));
           yakl::atomicAdd(crm_output_mcudn(l,icrm) , tmp);
         }
      }
      yakl::atomicAdd(crm_output_gliqwp(l,icrm) , qcl(k,j,i,icrm));
      yakl::atomicAdd(crm_output_gicewp(l,icrm) , qci(k,j,i,icrm));
    }
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // Reduced radiation method allows for fewer radiation calculations
    // by collecting statistics and doing radiation over column groups
    int i_rad = i / (nx/crm_nx_rad);
    int j_rad = j / (ny/crm_ny_rad);
    real qsat_tmp;
    real rh_tmp;

    yakl::atomicAdd(crm_rad_temperature(k,j_rad,i_rad,icrm) , tabs(k,j,i,icrm));
    real tmp = max(0.0,qv(k,j,i,icrm));
    yakl::atomicAdd(crm_rad_qv(k,j_rad,i_rad,icrm) , tmp);
    yakl::atomicAdd(crm_rad_qc(k,j_rad,i_rad,icrm) , qcl(k,j,i,icrm));
    yakl::atomicAdd(crm_rad_qi(k,j_rad,i_rad,icrm) , qci(k,j,i,icrm));
    if (qcl(k,j,i,icrm) + qci(k,j,i,icrm) > 0) {
      yakl::atomicAdd(crm_rad_cld(k,j_rad,i_rad,icrm) , CF3D(k,j,i,icrm));
    } else {
      qsatw_crm(tabs(k,j,i,icrm),pres(k,icrm),qsat_tmp);
      rh_tmp = qv(k,j,i,icrm)/qsat_tmp;
      yakl::atomicAdd(crm_clear_rh(k,icrm) , rh_tmp);
      yakl::atomicAdd(crm_clear_rh_cnt(k,icrm),1);
    }
  });

  // Diagnose mass fluxes to drive CAM's convective transport of tracers.
  // definition of mass fluxes is taken from Xu et al., 2002, QJRMS.

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  //      for (int k=0; k<nzm+1; k++) {
  parallel_for( SimpleBounds<4>(nzm+1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int l=plev+1-(k+1);
    int kx;
    real qsat;
    real tmp;
    if (w(k,j+offy_w,i+offx_w,icrm) > 0.0) {
      kx=max(0, k-1);
      qsatw_crm(tabs(kx,j,i,icrm),pres(kx,icrm),qsat);
      if (qcl(kx,j,i,icrm)+qci(kx,j,i,icrm) > min(1.0e-5,0.01*qsat)) {
        tmp = rhow(k,icrm)*w(k,j+offy_w,i+offx_w,icrm);
        yakl::atomicAdd(mui_crm(l,icrm) , tmp);
      }
    } else if (w(k,j+offy_w,i+offx_w,icrm) < 0.0) {
      kx=min(k+1, nzm-1);
      qsatw_crm(tabs(kx,j,i,icrm),pres(kx,icrm),qsat);
      if (qcl(kx,j,i,icrm)+qci(kx,j,i,icrm) > min(1.0e-5,0.01*qsat)) {
        tmp = rhow(k,icrm)*w(k,j+offy_w,i+offx_w,icrm);
        yakl::atomicAdd(mdi_crm(l,icrm) , tmp);
      } else if (qpl(kx,j,i,icrm)+qpi(kx,j,i,icrm) > 1.0e-4) {
        tmp = rhow(k,icrm)*w(k,j+offy_w,i+offx_w,icrm);
        yakl::atomicAdd(mdi_crm(l,icrm) , tmp);
      }
    }
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    if(cwp(j,i,icrm) > cwp_threshold) {
      yakl::atomicAdd(crm_output_cltot(icrm) , cttemp(j,i,icrm));
    }
    if(cwph(j,i,icrm) > cwp_threshold) {
      yakl::atomicAdd(crm_output_clhgh(icrm) , chtemp(j,i,icrm));
    }
    if(cwpm(j,i,icrm) > cwp_threshold) {
      yakl::atomicAdd(crm_output_clmed(icrm) , cmtemp(j,i,icrm));
    }
    if(cwpl(j,i,icrm) > cwp_threshold) {
      yakl::atomicAdd(crm_output_cllow(icrm) , cltemp(j,i,icrm));
    }
  });

}


