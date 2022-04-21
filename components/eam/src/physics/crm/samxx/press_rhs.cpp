#include "press_rhs.h"

void press_rhs() {
  YAKL_SCOPE( dudt          , :: dudt );
  YAKL_SCOPE( dvdt          , :: dvdt );
  YAKL_SCOPE( dwdt          , :: dwdt );
  YAKL_SCOPE( dx            , :: dx ); 
  YAKL_SCOPE( dy            , :: dy ); 
  YAKL_SCOPE( at            , :: at );
  YAKL_SCOPE( bt            , :: bt );
  YAKL_SCOPE( ct            , :: ct );
  YAKL_SCOPE( adz           , :: adz ); 
  YAKL_SCOPE( dz            , :: dz ); 
  YAKL_SCOPE( rhow          , :: rhow ); 
  YAKL_SCOPE( rho           , :: rho ); 
  YAKL_SCOPE( dt3           , :: dt3 );
  YAKL_SCOPE( na            , :: na ); 
  YAKL_SCOPE( nb            , :: nb );
  YAKL_SCOPE( nc            , :: nc );
  YAKL_SCOPE( u             , :: u );
  YAKL_SCOPE( v             , :: v );
  YAKL_SCOPE( w             , :: w );
  YAKL_SCOPE( p             , :: p );
  YAKL_SCOPE( ncrms         , :: ncrms );
  
  if (dowallx && rank%nsubdomains_x == 0) {
    // for (int k=0; k<nzm; k++) {
    //  for (int j=0; j<ny; j++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
      dudt(na-1,k,j,0,icrm) = 0.0;
    });
  }

  if (dowally && RUN3D == 1 && rank < nsubdomains_x) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dvdt(na-1,k,0,i,icrm) = 0.0;
    });
  }

  bound_duvdt();

  real rdx=1.0/dx;
  real rdy=1.0/dy;
  real btat=bt/at;
  real ctat=ct/at;

  if (RUN3D) {

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      real rdz=1.0/(adz(k,icrm)*dz(icrm));
      real rup = rhow(kc,icrm)/rho(k,icrm)*rdz;
      real rdn = rhow(k,icrm)/rho(k,icrm)*rdz;
      int jc=j+1;
      int ic=i+1;
      real dta=1.0/dt3(na-1)/at;
      p(k,j+offy_p,i+offx_p,icrm)=( rdx*(u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))+
                                  rdy*(v(k,jc+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))+
                                  (w(kc,j+offy_w,i+offx_w,icrm)*rup-w(k,j+offy_w,i+offx_w,icrm)*rdn) )*dta +
                                  ( rdx*(dudt(na-1,k,j,ic,icrm)-dudt(na-1,k,j,i,icrm))+
                                  rdy*(dvdt(na-1,k,jc,i,icrm)-dvdt(na-1,k,j,i,icrm))+
                                  (dwdt(na-1,kc,j,i,icrm)*rup-dwdt(na-1,k,j,i,icrm)*rdn) ) +
                                  btat*( rdx*(dudt(nb-1,k,j,ic,icrm)-dudt(nb-1,k,j,i,icrm))+
                                  rdy*(dvdt(nb-1,k,jc,i,icrm)-dvdt(nb-1,k,j,i,icrm))+
                                  (dwdt(nb-1,kc,j,i,icrm)*rup-dwdt(nb-1,k,j,i,icrm)*rdn) ) +
                                  ctat*( rdx*(dudt(nc-1,k,j,ic,icrm)-dudt(nc-1,k,j,i,icrm))+
                                  rdy*(dvdt(nc-1,k,jc,i,icrm)-dvdt(nc-1,k,j,i,icrm))+
                                  (dwdt(nc-1,kc,j,i,icrm)*rup-dwdt(nc-1,k,j,i,icrm)*rdn) );
      p(k,j+offy_p,i+offx_p,icrm)=p(k,j+offy_p,i+offx_p,icrm)*rho(k,icrm);
    });

  } else {

    int j=0;
    // for (int k=0; k<nzm; k++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=k+1;
      real rdz=1.0/(adz(k,icrm)*dz(icrm));
      real rup = rhow(kc,icrm)/rho(k,icrm)*rdz;
      real rdn = rhow(k,icrm)/rho(k,icrm)*rdz;
      int ic=i+1;
      real dta=1.0/dt3(na-1)/at;

      p(k,j+offy_p,i+offx_p,icrm)=(rdx*(u(k,j+offy_u,ic+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))+
                                  (w(kc,j+offy_w,i+offx_w,icrm)*rup-w(k,j+offy_w,i+offx_w,icrm)*rdn) )*dta +
                                  (rdx*(dudt(na-1,k,j,ic,icrm)-dudt(na-1,k,j,i,icrm))+
                                  (dwdt(na-1,kc,j,i,icrm)*rup-dwdt(na-1,k,j,i,icrm)*rdn) ) +
                                  btat*(rdx*(dudt(nb-1,k,j,ic,icrm)-dudt(nb-1,k,j,i,icrm))+
                                  (dwdt(nb-1,kc,j,i,icrm)*rup-dwdt(nb-1,k,j,i,icrm)*rdn) ) +
                                  ctat*(rdx*(dudt(nc-1,k,j,ic,icrm)-dudt(nc-1,k,j,i,icrm))+
                                  (dwdt(nc-1,kc,j,i,icrm)*rup-dwdt(nc-1,k,j,i,icrm)*rdn) );
                                  p(k,j+offy_p,i+offx_p,icrm)=p(k,j+offy_p,i+offx_p,icrm)*rho(k,icrm);
    });

  }
}


