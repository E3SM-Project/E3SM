#include "advect2_mom_z.h"

void advect2_mom_z() {

  auto &uwle          = :: uwle;
  auto &vwle          = :: vwle;
  auto &dz            = :: dz;
  auto &rhow          = :: rhow;
  auto &rho           = :: rho;
  auto &u             = :: u;
  auto &v             = :: v;
  auto &w             = :: w;
  auto &na            = :: na;
  auto &dudt          = :: dudt;
  auto &dvdt          = :: dvdt;
  auto &dwdt          = :: dwdt;
  auto &adz           = :: adz;
  auto &adzw          = :: adzw;
  auto &ncrms         = :: ncrms;

  real4d fuz("fuz",nz ,ny,nx,ncrms);
  real4d fvz("fvz",nz ,ny,nx,ncrms);
  real4d fwz("fwz",nzm,ny,nx,ncrms);

  // for (int k=0; k<nzm; k++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    uwle(k,icrm) = 0.0;
    vwle(k,icrm) = 0.0;
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    real dz25=1.0/(4.0*dz(icrm));
    fuz(0,j,i,icrm) = 0.0;
    fuz(nz-1,j,i,icrm) = 0.0;
    fvz(0,j,i,icrm) = 0.0;
    fvz(nz-1,j,i,icrm) = 0.0;
    fwz(0,j,i,icrm) = 0.0;
    fwz(nzm-1,j,i,icrm) = 0.0;
  });

  if (RUN3D) {

    // for (int k=0; k<nzm-1; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      real dz25=1.0/(4.0*dz(icrm));
      int kb = k-1;
      real rhoi = dz25 * rhow(k+1,icrm);
      fuz(k+1,j,i,icrm) = rhoi*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k+1,j+offy_w,i-1+offx_w,icrm))*
                               (u(k+1,j+offy_u,i+offx_u,icrm)+u(kb+1,j+offy_u,i+offx_u,icrm));
      fvz(k+1,j,i,icrm) = rhoi*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k+1,j-1+offy_w,i+offx_w,icrm))*
                               (v(k+1,j+offy_v,i+offx_v,icrm)+v(kb+1,j+offy_v,i+offx_v,icrm));
      yakl::atomicAdd(uwle(k+1,icrm),fuz(k+1,j,i,icrm));
      yakl::atomicAdd(vwle(k+1,icrm),fvz(k+1,j,i,icrm));
    });

  } else {

    // for (int k=0; k<nzm-1; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      real dz25=1.0/(4.0*dz(icrm));
      int kb = k-1;
      real rhoi = dz25 * rhow(k+1,icrm);
      real www = rhoi*(w(k+1,j+offy_w,i+offx_w,icrm)+w(k+1,j+offy_w,i-1+offx_w,icrm));
      fuz(k+1,j,i,icrm) = www*(u(k+1,j+offy_u,i+offx_u,icrm)+u(kb+1,j+offy_u,i+offx_u,icrm));
      fvz(k+1,j,i,icrm) = www*(v(k+1,j+offy_v,i+offx_v,icrm)+v(kb+1,j+offy_v,i+offx_v,icrm));
      yakl::atomicAdd(uwle(k+1,icrm),fuz(k+1,j,i,icrm));
      yakl::atomicAdd(vwle(k+1,icrm),fvz(k+1,j,i,icrm));
    });

  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real dz25=1.0/(4.0*dz(icrm));
    int kc = k+1;
    real rhoi = 1.0/(rho(k,icrm)*adz(k,icrm));
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fuz(kc,j,i,icrm)-fuz(k,j,i,icrm))*rhoi;
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fvz(kc,j,i,icrm)-fvz(k,j,i,icrm))*rhoi;
    fwz(k,j,i,icrm)=dz25*(w(kc,j+offy_w,i+offx_w,icrm)*rhow(kc,icrm)+w(k,j+offy_w,i+offx_w,icrm)*
                          rhow(k,icrm))*(w(kc,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,i+offx_w,icrm));
  });

  // for (int k=0; k<nzm-1; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kb=k-1;
    real rhoi = 1.0/(rhow(k+1,icrm)*adzw(k+1,icrm));
    dwdt(na-1,k+1,j,i,icrm)=dwdt(na-1,k+1,j,i,icrm)-(fwz(k+1,j,i,icrm)-fwz(kb+1,j,i,icrm))*rhoi;
  });

}
