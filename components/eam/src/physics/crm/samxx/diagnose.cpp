#include "diagnose.h"

void diagnose() {
  YAKL_SCOPE( u0             , ::u0);
  YAKL_SCOPE( v0             , ::v0);
  YAKL_SCOPE( t01            , ::t01);
  YAKL_SCOPE( q01            , ::q01);
  YAKL_SCOPE( t0             , ::t0);
  YAKL_SCOPE( tabs0          , ::tabs0);
  YAKL_SCOPE( q0             , ::q0);
  YAKL_SCOPE( qn0            , ::qn0);
  YAKL_SCOPE( qp0            , ::qp0);
  YAKL_SCOPE( p0             , ::p0);
  YAKL_SCOPE( rho            , ::rho);
  YAKL_SCOPE( dz             , ::dz);
  YAKL_SCOPE( adz            , ::adz);
  YAKL_SCOPE( dtfactor       , ::dtfactor);
  YAKL_SCOPE( tabs           , ::tabs);
  YAKL_SCOPE( t              , ::t);
  YAKL_SCOPE( gamaz          , ::gamaz); 
  YAKL_SCOPE( qcl            , ::qcl);
  YAKL_SCOPE( qpl            , ::qpl);
  YAKL_SCOPE( qci            , ::qci);
  YAKL_SCOPE( qpi            , ::qpi);
  YAKL_SCOPE( u              , ::u);
  YAKL_SCOPE( v              , ::v);
  YAKL_SCOPE( p              , ::p);
  YAKL_SCOPE( qv             , ::qv);
  YAKL_SCOPE( usfc_xy        , ::usfc_xy);
  YAKL_SCOPE( vsfc_xy        , ::vsfc_xy);
  YAKL_SCOPE( pres           , ::pres);
  YAKL_SCOPE( swvp_xy        , ::swvp_xy);
  YAKL_SCOPE( psfc_xy        , ::psfc_xy);
  YAKL_SCOPE( cloudtopheight , ::cloudtopheight);
  YAKL_SCOPE( cloudtoptemp   , ::cloudtoptemp);
  YAKL_SCOPE( echotopheight  , ::echotopheight);
  YAKL_SCOPE( sstxy          , ::sstxy);
  YAKL_SCOPE( z              , ::z);
  YAKL_SCOPE( cld_xy         , ::cld_xy);
  YAKL_SCOPE( qv0            , ::qv0);
  YAKL_SCOPE( ncrms          , ::ncrms);

  real coef = 1.0/( (real) nx * (real) ny );

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    u0   (k,icrm)=0.0;
    v0   (k,icrm)=0.0;
    t01  (k,icrm) = tabs0(k,icrm);
    q01  (k,icrm) = q0   (k,icrm);
    t0   (k,icrm)=0.0;
    tabs0(k,icrm)=0.0;
    q0   (k,icrm)=0.0;
    qn0  (k,icrm)=0.0;
    qp0  (k,icrm)=0.0;
    p0   (k,icrm)=0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor;
    tabs(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm)-gamaz(k,icrm)+ fac_cond *
                       (qcl(k,j,i,icrm)+qpl(k,j,i,icrm)) + fac_sub *(qci(k,j,i,icrm) + qpi(k,j,i,icrm));
    yakl::atomicAdd(u0(k,icrm),u(k,j+offy_u,i+offx_u,icrm));
    yakl::atomicAdd(v0(k,icrm),v(k,j+offy_v,i+offx_v,icrm));
    yakl::atomicAdd(p0(k,icrm),p(k,j+offy_p,i+offx_p,icrm));
    yakl::atomicAdd(t0(k,icrm),t(k,j+offy_s,i+offx_s,icrm));
    yakl::atomicAdd(tabs0(k,icrm),tabs(k,j,i,icrm));
    real tmp = qv(k,j,i,icrm)+qcl(k,j,i,icrm)+qci(k,j,i,icrm);
    yakl::atomicAdd(q0(k,icrm),tmp);
    tmp = qcl(k,j,i,icrm) + qci(k,j,i,icrm);
    yakl::atomicAdd(qn0(k,icrm),tmp);
    tmp = qpl(k,j,i,icrm) + qpi(k,j,i,icrm);
    yakl::atomicAdd(qp0(k,icrm),tmp);
    tmp = qv(k,j,i,icrm)*coef1;
    // TODO: There should seemingly be an atomic statement here
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    u0   (k,icrm)=u0   (k,icrm)*coef;
    v0   (k,icrm)=v0   (k,icrm)*coef;
    t0   (k,icrm)=t0   (k,icrm)*coef;
    tabs0(k,icrm)=tabs0(k,icrm)*coef;
    q0   (k,icrm)=q0   (k,icrm)*coef;
    qn0  (k,icrm)=qn0  (k,icrm)*coef;
    qp0  (k,icrm)=qp0  (k,icrm)*coef;
    p0   (k,icrm)=p0   (k,icrm)*coef;
  });

  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    usfc_xy(j,i,icrm) = usfc_xy(j,i,icrm) + u(0,j+offy_s,i+offx_s,icrm)*dtfactor;
    vsfc_xy(j,i,icrm) = vsfc_xy(j,i,icrm) + v(0,j+offy_s,i+offx_s,icrm)*dtfactor;
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qv0(k,icrm) = q0(k,icrm) - qn0(k,icrm);
  });

  //=====================================================
  // UW ADDITIONS
  // FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor;
    // Saturated water vapor path with respect to water. Can be used
    // with water vapor path (= pw) to compute column-average
    // relative humidity.
    real tmp;
    qsatw_crm(tabs(k,j,i,icrm),pres(k,icrm),tmp);
    tmp = tmp*coef1;
    yakl::atomicAdd(swvp_xy(j,i,icrm),tmp);
  });

  // ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS

  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    psfc_xy(j,i,icrm) = psfc_xy(j,i,icrm) + (100.0*pres(0,icrm) + p(0,j+offy_p,i+offx_p,icrm))*dtfactor;
  });

  // COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
  // WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
  // CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
  // WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.
  // initially, zero out heights and set cloudtoptemp to SST

  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    cloudtopheight(j,i,icrm) = 0.0;
    cloudtoptemp(j,i,icrm) = sstxy(j+offy_sstxy,i+offx_sstxy,icrm);
    echotopheight(j,i,icrm) = 0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    // FIND CLOUD TOP HEIGHT
    real tmp_lwp = 0.0;
    for(int k=nzm-1; k>=0; k--) {
      tmp_lwp = tmp_lwp + (qcl(k,j,i,icrm)+qci(k,j,i,icrm))*rho(k,icrm)*dz(icrm)*adz(k,icrm);
      if (tmp_lwp > 0.01) {
        cloudtopheight(j,i,icrm) = z(k,icrm);
        cloudtoptemp(j,i,icrm) = tabs(k,j,i,icrm);
        cld_xy(j,i,icrm) = cld_xy(j,i,icrm) + dtfactor;
        break;
      }
    }
    // FIND ECHO TOP HEIGHT
    for(int k=nzm-1; k>=0; k--) {
      if (qpl(k,j,i,icrm)+qpi(k,j,i,icrm) > 1.e-6) {
        echotopheight(j,i,icrm) = z(k,icrm);
        break;
      }
    }
  });
  // END UW ADDITIONS
  //=====================================================

}
