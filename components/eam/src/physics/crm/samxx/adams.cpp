
#include "adams.h"

void adams() {
  auto &dtn    = ::dtn   ;
  auto &dx     = ::dx    ;
  auto &dy     = ::dy    ;
  auto &dz     = ::dz    ;
  auto &rho    = ::rho   ;
  auto &rhow   = ::rhow  ;
  auto &dudt   = ::dudt  ;
  auto &dvdt   = ::dvdt  ;
  auto &dwdt   = ::dwdt  ;
  auto &u      = ::u     ;
  auto &v      = ::v     ;
  auto &w      = ::w     ;
  auto &misc   = ::misc  ;
  auto &dt3    = ::dt3   ;
  auto &na     = ::na    ;
  auto &nb     = ::nb    ;
  auto &nc     = ::nc    ;
  auto &at     = ::at    ;
  auto &bt     = ::bt    ;
  auto &ct     = ::ct    ;
  auto &ncrms  = ::ncrms ;

  // Adams-Bashforth scheme
  real dtdx = dtn/dx;
  real dtdy = dtn/dy;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real dtdz = dtn/dz(icrm);
    real rhox = rho (k,icrm)*dtdx;
    real rhoy = rho (k,icrm)*dtdy;
    real rhoz = rhow(k,icrm)*dtdz;
    real utend = ( at*dudt(na-1,k,j,i,icrm) + bt*dudt(nb-1,k,j,i,icrm) + ct*dudt(nc-1,k,j,i,icrm) );
    real vtend = ( at*dvdt(na-1,k,j,i,icrm) + bt*dvdt(nb-1,k,j,i,icrm) + ct*dvdt(nc-1,k,j,i,icrm) );
    real wtend = ( at*dwdt(na-1,k,j,i,icrm) + bt*dwdt(nb-1,k,j,i,icrm) + ct*dwdt(nc-1,k,j,i,icrm) );
    dudt(nc-1,k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm) + dt3(na-1) * utend;
    dvdt(nc-1,k,j,i,icrm) = v(k,j+offy_v,i+offx_v,icrm) + dt3(na-1) * vtend;
    dwdt(nc-1,k,j,i,icrm) = w(k,j+offy_w,i+offx_w,icrm) + dt3(na-1) * wtend;
    u   (k,j+offy_u,i+offx_u,icrm) = 0.5 * ( u(k,j+offy_u,i+offx_u,icrm) + dudt(nc-1,k,j,i,icrm) ) * rhox;
    v   (k,j+offy_v,i+offx_v,icrm) = 0.5 * ( v(k,j+offy_v,i+offx_v,icrm) + dvdt(nc-1,k,j,i,icrm) ) * rhoy;
    w   (k,j+offy_w,i+offx_w,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) ) * rhoz;
    misc(k,j       ,i       ,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) );
  });

}

