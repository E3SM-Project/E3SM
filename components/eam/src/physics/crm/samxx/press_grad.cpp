
#include "press_grad.h"

void press_grad() {
  YAKL_SCOPE( dz            , :: dz );
  YAKL_SCOPE( adzw          , :: adzw );
  YAKL_SCOPE( p             , :: p );
  YAKL_SCOPE( na            , :: na );
  YAKL_SCOPE( dudt          , :: dudt );
  YAKL_SCOPE( dvdt          , :: dvdt );
  YAKL_SCOPE( dwdt          , :: dwdt );
  YAKL_SCOPE( rho           , :: rho );
  YAKL_SCOPE( ncrms         , :: ncrms );

  real rdx=1.0/dx;
  real rdy=1.0/dy;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kb=max(0,k-1);
    real rdz = 1.0/(dz(icrm)*adzw(k,icrm));
    int jb=j-YES3D;
    int ib=i-1;
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(p(k,j+offy_p,i+offx_p,icrm)-p(k,j+offy_p,ib+offx_p,icrm))*rdx;
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(p(k,j+offy_p,i+offx_p,icrm)-p(k,jb+offy_p,i+offx_p,icrm))*rdy;
    dwdt(na-1,k,j,i,icrm)=dwdt(na-1,k,j,i,icrm)-(p(k,j+offy_p,i+offx_p,icrm)-p(kb,j+offy_p,i+offx_p,icrm))*rdz;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+YES3D,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    p(k,j,i,icrm)=p(k,j,i,icrm)*rho(k,icrm);  // convert p'/rho to p'
  });

  if(dowallx && rank%nsubdomains_x == 0) {
    // for (int k=0; k<nzm; k++) {
    //  for (int j=0; j<ny; j++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
      dudt(na-1,k,j,0,icrm) = 0.0;
    });
  }

  if(dowally && RUN3D == 1 && rank < nsubdomains_x) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dvdt(na-1,k,0,i,icrm) = 0.0;
    });
  }

  bound_duvdt();
}


