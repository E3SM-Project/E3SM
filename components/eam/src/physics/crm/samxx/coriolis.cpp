#include "coriolis.h"

void coriolis() {
  YAKL_SCOPE( u       , ::u);
  YAKL_SCOPE( v       , ::v);
  YAKL_SCOPE( w       , ::w);
  YAKL_SCOPE( fcory   , ::fcory);
  YAKL_SCOPE( fcorzy  , ::fcorzy);
  YAKL_SCOPE( dudt    , ::dudt);
  YAKL_SCOPE( dvdt    , ::dvdt);
  YAKL_SCOPE( dwdt    , ::dwdt);
  YAKL_SCOPE( na      , ::na);
  YAKL_SCOPE( vg0     , ::vg0);
  YAKL_SCOPE( ug0     , ::ug0);
  YAKL_SCOPE( ncrms   , ::ncrms);
  YAKL_SCOPE( use_ESMT, :: use_ESMT );
  YAKL_SCOPE( u_esmt  , :: u_esmt );
  YAKL_SCOPE( v_esmt  , :: v_esmt );

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
      // dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)+fcory(j+offy_fcory,icrm)*(v_av-vg0(k,icrm))-fcorzy(j,icrm)*w_av;
      // real u_av=0.25*(u(k,j+offy_u,i+offx_u,icrm)+u(k,j+offy_u,ic+offx_u,icrm)+u(k,jb+offy_u,i+offx_u,icrm)+
      //                 u(k,jb+offy_u,ic+offx_u,icrm));
      // dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-0.5*(fcory(j+offy_fcory,icrm)+fcory(jb+offy_fcory,icrm))*
      //                       (u_av-ug0(k,icrm));
      #ifdef MMF_DO_CORIOLIS_U
        // dudt(na-1,k,j,i,icrm) += fcory(j+offy_fcory,icrm)*(v_av-vg0(k,icrm));
        dudt(na-1,k,j,i,icrm) += -1*fcorzy(j,icrm) * w_av;
      #endif
      #ifdef MMF_DO_CORIOLIS_W             
        real u_av_on_w = 0.25*( u(k ,j+offy_u,i+offx_u,icrm) + u(k ,j+offy_u,ic+offx_u,icrm)
                               +u(kc,j+offy_u,i+offx_u,icrm) + u(kc,j+offy_u,ic+offx_u,icrm));
        dwdt(na-1,kc,j,i,icrm) += fcorzy(j,icrm)*u_av_on_w;
      #endif
    });
  } else {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    //   int kc=k+1;
    //   int kb=k-1;
    //   int ib=i-1;
    //   int ic=i+1;
      // real w_av=0.25*(w(kc,j+offy_w,i+offx_w,icrm)+w(kc,j+offy_w,ib+offx_w,icrm)+
      //                 w(k,j+offy_w,i+offx_w,icrm)+w(k,j+offy_w,ib+offx_w,icrm));
      // dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)+fcory(j+offy_fcory,icrm)*(v(k,j+offy_v,i+offx_v,icrm)-vg0(k,icrm))-fcorzy(j,icrm)*w_av;
      // dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-fcory(j+offy_fcory,icrm)*(u(k,j+offy_u,i+offx_u,icrm)-ug0(k,icrm));
    // });

    // #ifdef MMF_DO_CORIOLIS_W
    //   // calculate zonal wind perturbation
    //   real4d up_av_on_w("up_av_on_w",nzm,ny,nx,ncrms);
    //   real2d ub_av_on_w("ub_av_on_w",nzm,ncrms);
    //   real factor_xy = 1.0/( (real) nx * (real) ny );
      
      // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      //   int ic=i+1;
      //   int kc=k+1;
      //   real tmp;
      //   tmp = u(k ,j+offy_u,i+offx_u,icrm);
      // });
      // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      //   int ic=i+1;
      //   int kc=k+1;
      //   real tmp;
      //   tmp = u(k ,j+offy_u,ic+offx_u,icrm);
      // });
      // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      //   int ic=i+1;
      //   int kc=k+1;
      //   real tmp;
      //   tmp = u(kc,j+offy_u,i+offx_u,icrm);
      // });
      // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      //   int ic=i+1;
      //   int kc=k+1;
      //   real tmp;
      //   tmp = u(kc,j+offy_u,ic+offx_u,icrm);
      // });
      // parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      //   real tmp;
      //   tmp = up_av_on_w(k,j,i,icrm);
      // });

    #ifdef MMF_DO_CORIOLIS_W
      // calculate zonal wind perturbation
      real4d u_av_on_w("u_av_on_w",nzm+1,ny,nx,ncrms);
      real2d u_av_on_w_mean("u_av_on_w_mean",nzm+1,ncrms);
      
      parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        u_av_on_w(k+1,j,i,icrm) = 0.25*( u(k  ,j+offy_u,i  +offx_u,icrm)
                                        +u(k  ,j+offy_u,i+1+offx_u,icrm)
                                        +u(k+1,j+offy_u,i  +offx_u,icrm)
                                        +u(k+1,j+offy_u,i+1+offx_u,icrm));
      });
      
      real factor_xy = 1.0/( (real) nx * (real) ny );
      parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        yakl::atomicAdd( u_av_on_w_mean(k+1,icrm), u_av_on_w(k+1,j,i,icrm) / factor_xy );
      });

      parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        dwdt(na-1,k+1,j,i,icrm) += fcorzy(j,icrm) * ( u_av_on_w(k+1,j,i,icrm) - u_av_on_w_mean(k+1,icrm) );
      });
    #endif

    // ignore V wind terms for testing non-traditional coriolis terms
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=k+1;
      int kb=k-1;
      int ib=i-1;
      int ic=i+1;
      #ifdef MMF_DO_CORIOLIS_U
        real w_av_on_u = 0.25*( w(kc,j+offy_w,i+offx_w,icrm) + w(kc,j+offy_w,ib+offx_w,icrm)
                               +w(k ,j+offy_w,i+offx_w,icrm) + w(k ,j+offy_w,ib+offx_w,icrm));
        dudt(na-1,k,j,i,icrm) += -1*fcorzy(j,icrm)*w_av_on_u;
      #endif
      // #ifdef MMF_DO_CORIOLIS_W
        // if (kc<=nzm-1)
        // real u_av_on_w = 0.25*( u(k ,j+offy_u,i+offx_u,icrm) + u(k ,j+offy_u,ic+offx_u,icrm)
        //                        +u(kc,j+offy_u,i+offx_u,icrm) + u(kc,j+offy_u,ic+offx_u,icrm));
        // dwdt(na-1,kc,j,i,icrm) += fcorzy(j,icrm)*u_av_on_w;
        // dwdt(na-1,kc,j,i,icrm) += fcorzy(j,icrm)*up_av_on_w(k,j,i,icrm);
      // #endif
      #ifdef MMF_DO_CORIOLIS_ESMT
        if (use_ESMT) {
          real w_av_on_s = 0.5*( w(kc,j+offy_w,i+offx_w,icrm) + w(k,j+offy_w,i+offx_w,icrm) );
          u_esmt(k,j+offy_s,i+offx_s,icrm) += -1*fcorzy(j,icrm)*w_av_on_s;
        }
      #endif
      #ifdef MMF_DO_CORIOLIS_ESMT_W
        if (use_ESMT) {
          real u_esmt_av_on_w = 0.5*( u_esmt(kb,j+offy_s,i+offx_s,icrm) + u_esmt(k,j+offy_s,i+offx_s,icrm) );
          dwdt(na-1,kc,j,i,icrm) += fcorzy(j,icrm)*u_esmt_av_on_w;
        }
      #endif
    });
  }

}
