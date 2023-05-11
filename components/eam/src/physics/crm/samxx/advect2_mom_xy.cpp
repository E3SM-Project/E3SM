#include "advect2_mom_xy.h"

void advect2_mom_xy() {
  YAKL_SCOPE( dx             , :: dx);
  YAKL_SCOPE( dy             , :: dy);
  YAKL_SCOPE( rhow           , :: rhow); 
  YAKL_SCOPE( adzw           , :: adzw);
  YAKL_SCOPE( u              , :: u);
  YAKL_SCOPE( dudt           , :: dudt);
  YAKL_SCOPE( v              , :: v);
  YAKL_SCOPE( dvdt           , :: dvdt);
  YAKL_SCOPE( w              , :: w);
  YAKL_SCOPE( dwdt           , :: dwdt);
  YAKL_SCOPE( rho            , :: rho);
  YAKL_SCOPE( adz            , :: adz);
  YAKL_SCOPE( na             , :: na);
  YAKL_SCOPE( ncrms          , :: ncrms);

  real dx25 = 0.25 / dx;
  real dy25 = 0.25 / dy;

  if (RUN3D) {

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc= k+1;
      int kcu = min(kc, nzm-1);
      real irho = 1.0/(rhow(kc,icrm)*adzw(kc,icrm));
      int jb = j-1;
      int ic = i+1;
      real fu1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)+u(k,j+offy_u,i-1+offx_u,icrm))*
                      (u(k,j+offy_u,i-1+offx_u,icrm)+u(k,j+offy_u,ic-1+offx_u,icrm));
      real fu2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)+u(k,j+offy_u,i+offx_u,icrm))*
                      (u(k,j+offy_u,i+offx_u,icrm)+u(k,j+offy_u,ic+offx_u,icrm));
      dudt(na-1,k,j,i,icrm)  = dudt(na-1,k,j,i,icrm)  - (fu2-fu1);
      real fv1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)+u(k,jb+offy_u,ic-1+offx_u,icrm))*
                       (v(k,j+offy_v,i-1+offx_v,icrm)+v(k,j+offy_v,ic-1+offx_v,icrm));
      real fv2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)+u(k,jb+offy_u,ic+offx_u,icrm))*
                      (v(k,j+offy_v,i+offx_v,icrm)+v(k,j+offy_v,ic+offx_v,icrm));
      dvdt(na-1,k,j,i,icrm)  = dvdt(na-1,k,j,i,icrm)  - (fv2-fv1);
      real fw1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)*rho(k,icrm)*adz(k,icrm)+
                       u(kcu,j+offy_u,ic-1+offx_u,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*
                       (w(kc,j+offy_w,i-1+offx_w,icrm)+w(kc,j+offy_w,ic-1+offx_w,icrm));
      real fw2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)*rho(k,icrm)*adz(k,icrm)+
                       u(kcu,j+offy_u,ic+offx_u,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*
                       (w(kc,j+offy_w,i+offx_w,icrm)+w(kc,j+offy_w,ic+offx_w,icrm));
      dwdt(na-1,kc,j,i,icrm) = dwdt(na-1,kc,j,i,icrm)-irho*(fw2-fw1);

      int jc = j+1;
      int ib = i-1;
      fu1 = dy25*(v(k,jc-1+offy_v,i+offx_v,icrm)+v(k,jc-1+offy_v,ib+offx_v,icrm))*
                 (u(k,j-1+offy_u,i+offx_u,icrm)+u(k,jc-1+offy_u,i+offx_u,icrm));
      fu2 = dy25*(v(k,jc+offy_v,i+offx_v,icrm)+v(k,jc+offy_v,ib+offx_v,icrm))*
                 (u(k,j+offy_u,i+offx_u,icrm)+u(k,jc+offy_u,i+offx_u,icrm));
      dudt(na-1,k,j,i,icrm) = dudt(na-1,k,j,i,icrm) - (fu2-fu1);
      fv1 = dy25*(v(k,jc-1+offy_v,i+offx_v,icrm)+v(k,j-1+offy_v,i+offx_v,icrm))*
                 (v(k,j-1+offy_v,i+offx_v,icrm)+v(k,jc-1+offy_v,i+offx_v,icrm));
      fv2 = dy25*(v(k,jc+offy_v,i+offx_v,icrm)+v(k,j+offy_v,i+offx_v,icrm))*
                 (v(k,j+offy_v,i+offx_v,icrm)+v(k,jc+offy_v,i+offx_v,icrm));
      dvdt(na-1,k,j,i,icrm) = dvdt(na-1,k,j,i,icrm) - (fv2-fv1);
      fw1 = dy25*(v(k,jc-1+offy_v,i+offx_v,icrm)*rho(k,icrm)*adz(k,icrm)+
                  v(kcu,jc-1+offy_v,i+offx_v,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*
                  (w(kc,j-1+offy_w,i+offx_w,icrm)+w(kc,jc-1+offy_w,i+offx_w,icrm));
      fw2 = dy25*(v(k,jc+offy_v,i+offx_v,icrm)*rho(k,icrm)*adz(k,icrm)+
                  v(kcu,jc+offy_v,i+offx_v,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*
                  (w(kc,j+offy_w,i+offx_w,icrm)+w(kc,jc+offy_w,i+offx_w,icrm));
      dwdt(na-1,kc,j,i,icrm)= dwdt(na-1,kc,j,i,icrm)-irho*(fw2-fw1);
    });

  } else {

    // for (int k=0; k<nzm; k++) {
    //    for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int j=0;
      int kc= k+1;
      int kcu =min(kc, nzm-1);
      real irho = 1.0/(rhow(kc,icrm)*adzw(kc,icrm));
      int ic = i+1;
      real fu1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)+u(k,j+offy_u,i-1+offx_u,icrm))*
                      (u(k,j+offy_u,i-1+offx_u,icrm)+u(k,j+offy_u,ic-1+offx_u,icrm));
      real fu2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)+u(k,j+offy_u,i+offx_u,icrm))*
                      (u(k,j+offy_u,i+offx_u,icrm)+u(k,j+offy_u,ic+offx_u,icrm));
      dudt(na-1,k,j,i,icrm)  = dudt(na-1,k,j,i,icrm)  - (fu2-fu1);
      real fv1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)+u(k,j+offy_u,i-1+offx_u,icrm))*
                      (v(k,j+offy_v,i-1+offx_v,icrm)+v(k,j+offy_v,ic-1+offx_v,icrm));
      real fv2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)+u(k,j+offy_u,i+offx_u,icrm))*
                      (v(k,j+offy_v,i+offx_v,icrm)+v(k,j+offy_v,ic+offx_v,icrm));
      dvdt(na-1,k,j,i,icrm)  = dvdt(na-1,k,j,i,icrm)  - (fv2-fv1);
      real fw1 = dx25*(u(k,j+offy_u,ic-1+offx_u,icrm)*rho(k,icrm)*adz(k,icrm)+u(kcu,j+offy_u,ic-1+offx_u,icrm)*
                       rho(kcu,icrm)*adz(kcu,icrm))*(w(kc,j+offy_w,i-1+offx_w,icrm)+w(kc,j+offy_w,ic-1+offx_w,icrm));
      real fw2 = dx25*(u(k,j+offy_u,ic+offx_u,icrm)*rho(k,icrm)*adz(k,icrm)+u(kcu,j+offy_u,ic+offx_u,icrm)*
                       rho(kcu,icrm)*adz(kcu,icrm))*(w(kc,j+offy_w,i+offx_w,icrm)+w(kc,j+offy_w,ic+offx_w,icrm));
      dwdt(na-1,kc,j,i,icrm) = dwdt(na-1,kc,j,i,icrm)-irho*(fw2-fw1);
    });

  }

}
