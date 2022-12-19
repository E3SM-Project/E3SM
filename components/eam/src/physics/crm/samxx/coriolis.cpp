#include "coriolis.h"

void coriolis() {
  YAKL_SCOPE( u       , ::u);
  YAKL_SCOPE( v       , ::v);
  YAKL_SCOPE( w       , ::w);
  YAKL_SCOPE( fcory   , ::fcory);
  YAKL_SCOPE( fcorzy  , ::fcorzy);
  YAKL_SCOPE( dudt    , ::dudt);
  YAKL_SCOPE( dvdt    , ::dvdt);
  YAKL_SCOPE( na      , ::na);
  YAKL_SCOPE( vg0     , ::vg0);
  YAKL_SCOPE( ug0     , ::ug0);
  YAKL_SCOPE( ncrms   , ::ncrms);

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      real v_av=0.25*(v(k,j+offy_v,i+offx_v,icrm)+v(k,jc+offy_v,i+offx_v,icrm)+
                      v(k,j+offy_v,ib+offx_v,icrm)+v(k,jc+offy_v,ib+offx_v,icrm));
      real w_av=0.25*(w(kc,j+offy_w,i+offx_w,icrm)+w(kc,j+offy_w,ib+offx_w,icrm)+
                      w(k,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,ib+offx_w,icrm));
      dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)+fcory(j+offy_fcory,icrm)*(v_av-vg0(k,icrm))-fcorzy(j,icrm)*w_av;
      real u_av=0.25*(u(k,j+offy_u,i+offx_u,icrm)+u(k,j+offy_u,ic+offx_u,icrm)+u(k,jb+offy_u,i+offx_u,icrm)+
                      u(k,jb+offy_u,ic+offx_u,icrm));
      dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-0.5*(fcory(j+offy_fcory,icrm)+fcory(jb+offy_fcory,icrm))*
                            (u_av-ug0(k,icrm));
    });
  } else {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      int ib=i-1;
      int ic=i+1;
      real w_av=0.25*(w(kc,j+offy_w,i+offx_w,icrm)+w(kc,j+offy_w,ib+offx_w,icrm)+
                      w(k,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,ib+offx_w,icrm));
      dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)+fcory(j+offy_fcory,icrm)*
                            (v(k,j+offy_v,i+offx_v,icrm)-vg0(k,icrm))-fcorzy(j,icrm)*w_av;
      dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-fcory(j+offy_fcory,icrm)*
                            (u(k,j+offy_u,i+offx_u,icrm)-ug0(k,icrm));
    });
  }

}
