#include "diagnose.h"

void diagnose() {
  auto &u0             = ::u0;
  auto &v0             = ::v0;
  auto &t01            = ::t01;
  auto &q01            = ::q01;
  auto &t0             = ::t0;
  auto &tabs0          = ::tabs0;
  auto &q0             = ::q0;
  auto &qn0            = ::qn0;
  auto &qp0            = ::qp0;
  auto &p0             = ::p0;
  auto &rho            = ::rho;
  auto &dz             = ::dz;
  auto &adz            = ::adz;
  auto &dtfactor       = ::dtfactor;
  auto &tabs           = ::tabs;
  auto &t              = ::t;
  auto &gamaz          = ::gamaz; 
  auto &qcl            = ::qcl;
  auto &qpl            = ::qpl;
  auto &qci            = ::qci;
  auto &qpi            = ::qpi;
  auto &u              = ::u;
  auto &v              = ::v;
  auto &p              = ::p;
  auto &qv             = ::qv;
  auto &usfc_xy        = ::usfc_xy;
  auto &vsfc_xy        = ::vsfc_xy;
  auto &pres           = ::pres;
  auto &swvp_xy        = ::swvp_xy;
  auto &psfc_xy        = ::psfc_xy;
  auto &cloudtopheight = ::cloudtopheight;
  auto &cloudtoptemp   = ::cloudtoptemp;
  auto &echotopheight  = ::echotopheight;
  auto &sstxy          = ::sstxy;
  auto &z              = ::z;
  auto &cld_xy         = ::cld_xy;
  auto &qv0            = ::qv0;
  auto &ncrms          = ::ncrms;

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
