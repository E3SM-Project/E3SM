#include "advect_scalar3D.h"

void advect_scalar3D(real4d &f, real2d &flux) {
  YAKL_SCOPE( dowallx  , ::dowallx);
  YAKL_SCOPE( dowally  , ::dowally);
  YAKL_SCOPE( rank     , ::rank);
  YAKL_SCOPE( u        , ::u);
  YAKL_SCOPE( v        , ::v);
  YAKL_SCOPE( w        , ::w);
  YAKL_SCOPE( rho      , ::rho);
  YAKL_SCOPE( adz      , ::adz);
  YAKL_SCOPE( rhow     , ::rhow);
  YAKL_SCOPE( ncrms    , ::ncrms);

  bool constexpr nonos    = true;
  real constexpr eps      = 1.0e-10;
  int  constexpr offx_m   = 1;
  int  constexpr offy_m   = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offy_uuu = 2;
  int  constexpr offx_vvv = 2;
  int  constexpr offy_vvv = 2;
  int  constexpr offx_www = 2;
  int  constexpr offy_www = 2;

  real4d mx   ("mx"   ,nzm,ny+2,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,ny+2,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,ny+4,nx+5,ncrms);
  real4d vvv  ("vvv"  ,nzm,ny+5,nx+4,ncrms);
  real4d www  ("www"  ,nz ,ny+4,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<1-dimx1_u+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,1-dimx1_u+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x == nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,dimx2_u-(nx+1)+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int iInd = i+(nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (dowally) {
    if (rank < nsubdomains_x) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<1-dimy1_v+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,1-dimy1_v+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        v(k,j,i,icrm) = 0.0;
      });
    }
    if (rank > nsubdomains-nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy2_v-(ny+1)+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy2_v-(ny+1)+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int jInd = j+(ny+2);
        v(k,jInd,i,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
           max(f(k,j+offy_s-1,ib+offx_s-1,icrm),max(f(k,j+offy_s-1,ic+offx_s-1,icrm),
           max(f(k,jb+offy_s-1,i+offx_s-1,icrm),max(f(k,jc+offy_s-1,i+offx_s-1,icrm),
           max(f(kb,j+offy_s-1,i+offx_s-1,icrm),max(f(kc,j+offy_s-1,i+offx_s-1,icrm),f(k,j+offy_s-1,i+offx_s-1,icrm)))))));
      mn(k,j,i,icrm) = 
           min(f(k,j+offy_s-1,ib+offx_s-1,icrm),min(f(k,j+offy_s-1,ic+offx_s-1,icrm),
           min(f(k,jb+offy_s-1,i+offx_s-1,icrm),min(f(k,jc+offy_s-1,i+offx_s-1,icrm),
           min(f(kb,j+offy_s-1,i+offx_s-1,icrm),min(f(kc,j+offy_s-1,i+offx_s-1,icrm),f(k,j+offy_s-1,i+offx_s-1,icrm)))))));
    });
  } 

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+5; j++) {
  //     for (int i=0; i<nx+5; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+5,nx+5,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kb=max(0,k-1);
    if (j <= ny+3){
      uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(k,j+offy_s-2,i-1+offx_s-2,icrm)+
                      min(0.0,u(k,j,i,icrm))*f(k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i <= nx+3) {
      vvv(k,j,i,icrm)=max(0.0,v(k,j,i,icrm))*f(k,j-1+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,v(k,j,i,icrm))*f(k,j+offx_s-2,i+offy_s-2,icrm);
    }
    if (i <= nx+3 && j <= ny+3) {
      www(k,j,i,icrm)=max(0.0,w(k,j,i,icrm))*f(kb,j+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,w(k,j,i,icrm))*f(k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i == 0 && j == 0) {
      flux(k,icrm) = 0.0;
    }
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    irho(k,icrm) = 1.0/rho(k,icrm);
    iadz(k,icrm) = 1.0/adz(k,icrm);
    irhow(k,icrm) = 1.0/(rhow(k,icrm)*adz(k,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (i >= 2 && i <= nx+1 && j >= 2 && j <= ny+1) {
      yakl::atomicAdd(flux(k,icrm),www(k,j,i,icrm));
    }
    f(k,j+offy_s-2,i+offy_s-2,icrm)=f(k,j+offy_s-2,i+offx_s-2,icrm)-( uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                    vvv(k,j+1,i,icrm)-vvv(k,j,i,icrm)
                                    +(www(k+1,j,i,icrm)-www(k,j,i,icrm) )*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+3; j++) {
  //     for (int i=0; i<nx+3; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+3,nx+3,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (j <= ny+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm) = 
           andiff(f(k,j+offy_s-1,ib+offx_s-1,icrm),f(k,j+offy_s-1,i+offx_s-1,icrm),u(k,j+offy_u-1,i+offx_u-1,icrm),irho(k,icrm))-
          (across(f(k,jc+offy_s-1,ib+offx_s-1,icrm)+f(k,jc+offy_s-1,i+offx_s-1,icrm)-f(k,jb+offy_s-1,ib+offx_s-1,icrm)-
                  f(k,jb+offy_s-1,i+offx_s-1,icrm),u(k,j+offy_u-1,i+offx_u-1,icrm), v(k,j+offy_v-1,ib+offx_v-1,icrm)+
                  v(k,jc+offy_v-1,ib+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+v(k,j+offy_v-1,i+offx_v-1,icrm))+
           across(dd*(f(kc,j+offy_s-1,ib+offx_s-1,icrm)+f(kc,j+offy_s-1,i+offx_s-1,icrm)-f(kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(kb,j+offy_s-1,i+offx_s-1,icrm)),u(k,j+offy_u-1,i+offx_u-1,icrm), w(k,j+offy_w-1,ib+offx_w-1,icrm)+
                  w(kc,j+offy_w-1,ib+offx_w-1,icrm)+w(k,j+offy_w-1,i+offx_w-1,icrm)+w(kc,j+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int ib=i-1;
      int ic=i+1;
      vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm) = 
           andiff(f(k,jb+offy_s-1,i+offx_s-1,icrm),f(k,j+offy_s-1,i+offx_s-1,icrm),v(k,j+offy_v-1,i+offx_v-1,icrm),irho(k,icrm))-
           (across(f(k,jb+offy_s-1,ic+offx_s-1,icrm)+f(k,j+offy_s-1,ic+offx_s-1,icrm)-f(k,jb+offy_s-1,ib+offx_s-1,icrm)-
                   f(k,j+offy_s-1,ib+offx_s-1,icrm),v(k,j+offy_v-1,i+offx_v-1,icrm), u(k,jb+offy_u-1,i+offx_u-1,icrm)+
                   u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+u(k,jb+offy_u-1,ic+offx_u-1,icrm))+
            across(dd*(f(kc,jb+offy_s-1,i+offx_s-1,icrm)+f(kc,j+offy_s-1,i+offx_s-1,icrm)-f(kb,jb+offy_s-1,i+offx_s-1,icrm)-
                   f(kb,j+offy_s-1,i+offx_s-1,icrm)),v(k,j+offy_v-1,i+offx_v-1,icrm), w(k,jb+offy_w-1,i+offx_w-1,icrm)+
                   w(k,j+offy_w-1,i+offx_w-1,icrm)+w(kc,j+offy_w-1,i+offx_w-1,icrm)+w(kc,jb+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1 && j <= ny+1) {
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      www(k,j+offy_www-1,i+offx_www-1,icrm) = 
           andiff(f(kb,j+offy_s-1,i+offx_s-1,icrm),f(k,j+offy_s-1,i+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm),irhow(k,icrm))-
          (across(f(kb,j+offy_s-1,ic+offx_s-1,icrm)+f(k,j+offy_s-1,ic+offx_s-1,icrm)-f(kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(k,j+offy_s-1,ib+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), u(kb,j+offy_u-1,i+offx_u-1,icrm)+
                  u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+u(kb,j+offy_u-1,ic+offx_u-1,icrm))+
           across(f(k,jc+offy_s-1,i+offx_s-1,icrm)+f(kb,jc+offy_s-1,i+offx_s-1,icrm)-f(k,jb+offy_s-1,i+offx_s-1,icrm)-
                  f(kb,jb+offy_s-1,i+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), v(kb,j+offy_v-1,i+offx_v-1,icrm)+
                  v(kb,jc+offy_v-1,i+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+v(k,j+offy_v-1,i+offx_v-1,icrm))) *irho(k,icrm);
    }
  });

  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
          max(f(k,j+offy_s-1,ib+offx_s-1,icrm),max(f(k,j+offy_s-1,ic+offx_s-1,icrm),max(f(k,jb+offy_s-1,i+offx_s-1,icrm),
          max(f(k,jc+offy_s-1,i+offx_s-1,icrm),max(f(kb,j+offy_s-1,i+offx_s-1,icrm),max(f(kc,j+offy_s-1,i+offx_s-1,icrm),
          max(f(k,j+offy_s-1,i+offx_s-1,icrm),mx(k,j,i,icrm))))))));
      mn(k,j,i,icrm) = 
          min(f(k,j+offy_s-1,ib+offx_s-1,icrm),min(f(k,j+offy_s-1,ic+offx_s-1,icrm),min(f(k,jb+offy_s-1,i+offx_s-1,icrm),
          min(f(k,jc+offy_s-1,i+offx_s-1,icrm),min(f(kb,j+offy_s-1,i+offx_s-1,icrm),min(f(kc,j+offy_s-1,i+offx_s-1,icrm),
          min(f(k,j+offy_s-1,i+offx_s-1,icrm),mn(k,j,i,icrm))))))));
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int jc=j+1;
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(k,j+offy_s-1,i+offx_s-1,icrm))/
                ( pn3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pp3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pn3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pp3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pn3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pp3(www(k,j+offy_www-1,i+offx_www-1,icrm)))*iadz(k,icrm)+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(k,j+offy_s-1,i+offx_s-1,icrm)-mn(k,j,i,icrm))/
                ( pp3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pn3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pp3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pn3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pp3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pn3(www(k,j+offy_www-1,i+offx_www-1,icrm)))*iadz(k,icrm)+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j <= ny-1) {
        int ib=i-1;
        uuu(k,j+offy_uuu,i+offx_uuu,icrm) = 
              pp3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(k,j+offy_m,ib+offx_m,icrm)))
             -pn3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,ib+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1) {
        int jb=j-1;
        vvv(k,j+offy_vvv,i+offx_vvv,icrm) =
              pp3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(k,jb+offy_m,i+offx_m,icrm)))
             -pn3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,jb+offy_m,i+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1 && j <= ny-1) {
        int kb=max(0,k-1);
        www(k,j+offy_www,i+offx_www,icrm) =
              pp3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(kb,j+offy_m,i+offx_m,icrm)))
             -pn3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(kb,j+offy_m,i+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
        yakl::atomicAdd(flux(k,icrm),www(k,j+offy_www,i+offx_www,icrm));
      }
    });
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    int kc=k+1;
    f(k,j+offy_s,i+offx_s,icrm) = 
         max(0.0,f(k,j+offy_s,i+offx_s,icrm) -(uuu(k,j+offy_uuu,i+offx_uuu+1,icrm)-uuu(k,j+offy_uuu,i+offx_uuu,icrm)+
                 vvv(k,j+offy_vvv+1,i+offx_vvv,icrm)-vvv(k,j+offy_vvv,i+offx_vvv,icrm)+(www(k+1,j+offy_www,i+offx_www,icrm)-
                 www(k,j+offy_www,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}

void advect_scalar3D(real5d &f, int ind_f, real2d &flux) {
  YAKL_SCOPE( dowallx  , ::dowallx);
  YAKL_SCOPE( dowally  , ::dowally);
  YAKL_SCOPE( rank     , ::rank);
  YAKL_SCOPE( u        , ::u);
  YAKL_SCOPE( v        , ::v);
  YAKL_SCOPE( w        , ::w);
  YAKL_SCOPE( rho      , ::rho);
  YAKL_SCOPE( adz      , ::adz);
  YAKL_SCOPE( rhow     , ::rhow);
  YAKL_SCOPE( ncrms    , ::ncrms);

  bool constexpr nonos    = true;
  real constexpr eps      = 1.0e-10;
  int  constexpr offx_m   = 1;
  int  constexpr offy_m   = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offy_uuu = 2;
  int  constexpr offx_vvv = 2;
  int  constexpr offy_vvv = 2;
  int  constexpr offx_www = 2;
  int  constexpr offy_www = 2;

  real4d mx   ("mx"   ,nzm,ny+2,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,ny+2,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,ny+4,nx+5,ncrms);
  real4d vvv  ("vvv"  ,nzm,ny+5,nx+4,ncrms);
  real4d www  ("www"  ,nz ,ny+4,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<1-dimx1_u+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,1-dimx1_u+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x == nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,dimx2_u-(nx+1)+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int iInd = i+(nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (dowally) {
    if (rank < nsubdomains_x) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<1-dimy1_v+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,1-dimy1_v+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        v(k,j,i,icrm) = 0.0;
      });
    }
    if (rank > nsubdomains-nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy2_v-(ny+1)+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy2_v-(ny+1)+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int jInd = j+(ny+2);
        v(k,jInd,i,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
           max(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),max(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
           max(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),
           max(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)))))));
      mn(k,j,i,icrm) = 
           min(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),min(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
           min(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),
           min(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)))))));
    });
  } 

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+5; j++) {
  //     for (int i=0; i<nx+5; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+5,nx+5,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kb=max(0,k-1);
    if (j <= ny+3){
      uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i-1+offx_s-2,icrm)+
                      min(0.0,u(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i <= nx+3) {
      vvv(k,j,i,icrm)=max(0.0,v(k,j,i,icrm))*f(ind_f,k,j-1+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,v(k,j,i,icrm))*f(ind_f,k,j+offx_s-2,i+offy_s-2,icrm);
    }
    if (i <= nx+3 && j <= ny+3) {
      www(k,j,i,icrm)=max(0.0,w(k,j,i,icrm))*f(ind_f,kb,j+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,w(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i == 0 && j == 0) {
      flux(k,icrm) = 0.0;
    }
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    irho(k,icrm) = 1.0/rho(k,icrm);
    iadz(k,icrm) = 1.0/adz(k,icrm);
    irhow(k,icrm) = 1.0/(rhow(k,icrm)*adz(k,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (i >= 2 && i <= nx+1 && j >= 2 && j <= ny+1) {
      yakl::atomicAdd(flux(k,icrm),www(k,j,i,icrm));
    }
    f(ind_f,k,j+offy_s-2,i+offy_s-2,icrm)=f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm)-( uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                    vvv(k,j+1,i,icrm)-vvv(k,j,i,icrm)
                                    +(www(k+1,j,i,icrm)-www(k,j,i,icrm) )*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+3; j++) {
  //     for (int i=0; i<nx+3; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+3,nx+3,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (j <= ny+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm) = 
           andiff(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),u(k,j+offy_u-1,i+offx_u-1,icrm),irho(k,icrm))-
          (across(f(ind_f,k,jc+offy_s-1,ib+offx_s-1,icrm)+f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm)-f(ind_f,k,jb+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),u(k,j+offy_u-1,i+offx_u-1,icrm), v(k,j+offy_v-1,ib+offx_v-1,icrm)+
                  v(k,jc+offy_v-1,ib+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+v(k,j+offy_v-1,i+offx_v-1,icrm))+
           across(dd*(f(ind_f,kc,j+offy_s-1,ib+offx_s-1,icrm)+f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm)-f(ind_f,kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm)),u(k,j+offy_u-1,i+offx_u-1,icrm), w(k,j+offy_w-1,ib+offx_w-1,icrm)+
                  w(kc,j+offy_w-1,ib+offx_w-1,icrm)+w(k,j+offy_w-1,i+offx_w-1,icrm)+w(kc,j+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int ib=i-1;
      int ic=i+1;
      vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm) = 
           andiff(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),v(k,j+offy_v-1,i+offx_v-1,icrm),irho(k,icrm))-
           (across(f(ind_f,k,jb+offy_s-1,ic+offx_s-1,icrm)+f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm)-f(ind_f,k,jb+offy_s-1,ib+offx_s-1,icrm)-
                   f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),v(k,j+offy_v-1,i+offx_v-1,icrm), u(k,jb+offy_u-1,i+offx_u-1,icrm)+
                   u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+u(k,jb+offy_u-1,ic+offx_u-1,icrm))+
            across(dd*(f(ind_f,kc,jb+offy_s-1,i+offx_s-1,icrm)+f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm)-f(ind_f,kb,jb+offy_s-1,i+offx_s-1,icrm)-
                   f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm)),v(k,j+offy_v-1,i+offx_v-1,icrm), w(k,jb+offy_w-1,i+offx_w-1,icrm)+
                   w(k,j+offy_w-1,i+offx_w-1,icrm)+w(kc,j+offy_w-1,i+offx_w-1,icrm)+w(kc,jb+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1 && j <= ny+1) {
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      www(k,j+offy_www-1,i+offx_www-1,icrm) = 
           andiff(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm),irhow(k,icrm))-
          (across(f(ind_f,kb,j+offy_s-1,ic+offx_s-1,icrm)+f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm)-f(ind_f,kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), u(kb,j+offy_u-1,i+offx_u-1,icrm)+
                  u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+u(kb,j+offy_u-1,ic+offx_u-1,icrm))+
           across(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm)+f(ind_f,kb,jc+offy_s-1,i+offx_s-1,icrm)-f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm)-
                  f(ind_f,kb,jb+offy_s-1,i+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), v(kb,j+offy_v-1,i+offx_v-1,icrm)+
                  v(kb,jc+offy_v-1,i+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+v(k,j+offy_v-1,i+offx_v-1,icrm))) *irho(k,icrm);
    }
  });

  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
          max(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),max(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),max(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),
          max(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
          max(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),mx(k,j,i,icrm))))))));
      mn(k,j,i,icrm) = 
          min(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),min(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),min(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),
          min(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
          min(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),mn(k,j,i,icrm))))))));
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int jc=j+1;
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm))/
                ( pn3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pp3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pn3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pp3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pn3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pp3(www(k,j+offy_www-1,i+offx_www-1,icrm)))*iadz(k,icrm)+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)-mn(k,j,i,icrm))/
                ( pp3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pn3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pp3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pn3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pp3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pn3(www(k,j+offy_www-1,i+offx_www-1,icrm)))*iadz(k,icrm)+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j <= ny-1) {
        int ib=i-1;
        uuu(k,j+offy_uuu,i+offx_uuu,icrm) = 
              pp3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(k,j+offy_m,ib+offx_m,icrm)))
             -pn3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,ib+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1) {
        int jb=j-1;
        vvv(k,j+offy_vvv,i+offx_vvv,icrm) =
              pp3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(k,jb+offy_m,i+offx_m,icrm)))
             -pn3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,jb+offy_m,i+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1 && j <= ny-1) {
        int kb=max(0,k-1);
        www(k,j+offy_www,i+offx_www,icrm) =
              pp3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), mn(kb,j+offy_m,i+offx_m,icrm)))
             -pn3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(kb,j+offy_m,i+offx_m,icrm),mn(k,j+offy_m,i+offx_m,icrm)));
        yakl::atomicAdd(flux(k,icrm),www(k,j+offy_www,i+offx_www,icrm));
      }
    });
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    int kc=k+1;
    f(ind_f,k,j+offy_s,i+offx_s,icrm) = 
         max(0.0,f(ind_f,k,j+offy_s,i+offx_s,icrm) -(uuu(k,j+offy_uuu,i+offx_uuu+1,icrm)-uuu(k,j+offy_uuu,i+offx_uuu,icrm)+
                 vvv(k,j+offy_vvv+1,i+offx_vvv,icrm)-vvv(k,j+offy_vvv,i+offx_vvv,icrm)+(www(k+1,j+offy_www,i+offx_www,icrm)-
                 www(k,j+offy_www,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}

void advect_scalar3D(real5d &f, int ind_f, real3d &flux, int ind_flux) {
  YAKL_SCOPE( dowallx  , ::dowallx);
  YAKL_SCOPE( dowally  , ::dowally);
  YAKL_SCOPE( rank     , ::rank);
  YAKL_SCOPE( u        , ::u);
  YAKL_SCOPE( v        , ::v);
  YAKL_SCOPE( w        , ::w);
  YAKL_SCOPE( rho      , ::rho);
  YAKL_SCOPE( adz      , ::adz);
  YAKL_SCOPE( rhow     , ::rhow);
  YAKL_SCOPE( ncrms    , ::ncrms);

  bool constexpr nonos    = true;
  real constexpr eps      = 1.0e-10;
  int  constexpr offx_m   = 1;
  int  constexpr offy_m   = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offy_uuu = 2;
  int  constexpr offx_vvv = 2;
  int  constexpr offy_vvv = 2;
  int  constexpr offx_www = 2;
  int  constexpr offy_www = 2;

  real4d mx   ("mx"   ,nzm,ny+2,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,ny+2,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,ny+4,nx+5,ncrms);
  real4d vvv  ("vvv"  ,nzm,ny+5,nx+4,ncrms);
  real4d www  ("www"  ,nz ,ny+4,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<1-dimx1_u+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,1-dimx1_u+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x == nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy_u; j++) {
      //     for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy_u,dimx2_u-(nx+1)+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int iInd = i+(nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (dowally) {
    if (rank < nsubdomains_x) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<1-dimy1_v+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,1-dimy1_v+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        v(k,j,i,icrm) = 0.0;
      });
    }
    if (rank > nsubdomains-nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //   for (int j=0; j<dimy2_v-(ny+1)+1; j++) {
      //     for (int i=0; i<dimx_v; i++) {
      //       for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,dimy2_v-(ny+1)+1,dimx_v,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        int jInd = j+(ny+2);
        v(k,jInd,i,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
           max(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),max(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
           max(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),
           max(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
                                                          f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)))))));
      mn(k,j,i,icrm) = 
           min(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),min(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
           min(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),
           min(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
                                                          f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)))))));
    });
  } 

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+5; j++) {
  //     for (int i=0; i<nx+5; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+5,nx+5,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kb=max(0,k-1);
    if (j <= ny+3){
      uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i-1+offx_s-2,icrm)+
                      min(0.0,u(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i <= nx+3) {
      vvv(k,j,i,icrm)=max(0.0,v(k,j,i,icrm))*f(ind_f,k,j-1+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,v(k,j,i,icrm))*f(ind_f,k,j+offx_s-2,i+offy_s-2,icrm);
    }
    if (i <= nx+3 && j <= ny+3) {
      www(k,j,i,icrm)=max(0.0,w(k,j,i,icrm))*f(ind_f,kb,j+offy_s-2,i+offx_s-2,icrm)+
                      min(0.0,w(k,j,i,icrm))*f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm);
    }
    if (i == 0 && j == 0) {
      flux(ind_flux,k,icrm) = 0.0;
    }
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    irho(k,icrm) = 1.0/rho(k,icrm);
    iadz(k,icrm) = 1.0/adz(k,icrm);
    irhow(k,icrm) = 1.0/(rhow(k,icrm)*adz(k,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (i >= 2 && i <= nx+1 && j >= 2 && j <= ny+1) {
      yakl::atomicAdd(flux(ind_flux,k,icrm),www(k,j,i,icrm));
    }
    f(ind_f,k,j+offy_s-2,i+offy_s-2,icrm)=f(ind_f,k,j+offy_s-2,i+offx_s-2,icrm)-( uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                    vvv(k,j+1,i,icrm)-vvv(k,j,i,icrm)
                                    +(www(k+1,j,i,icrm)-www(k,j,i,icrm) )*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+3; j++) {
  //     for (int i=0; i<nx+3; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+3,nx+3,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (j <= ny+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm) = 
           andiff(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),
                  u(k,j+offy_u-1,i+offx_u-1,icrm),irho(k,icrm))-
          (across(f(ind_f,k,jc+offy_s-1,ib+offx_s-1,icrm)+f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm)-
                  f(ind_f,k,jb+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),u(k,j+offy_u-1,i+offx_u-1,icrm),
                  v(k,j+offy_v-1,ib+offx_v-1,icrm)+
                  v(k,jc+offy_v-1,ib+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+
                  v(k,j+offy_v-1,i+offx_v-1,icrm))+
           across(dd*(f(ind_f,kc,j+offy_s-1,ib+offx_s-1,icrm)+f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm)-
                  f(ind_f,kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm)),u(k,j+offy_u-1,i+offx_u-1,icrm), 
                  w(k,j+offy_w-1,ib+offx_w-1,icrm)+
                  w(kc,j+offy_w-1,ib+offx_w-1,icrm)+w(k,j+offy_w-1,i+offx_w-1,icrm)+
                  w(kc,j+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      real dd=2.0/(kc-kb)/adz(k,icrm);
      int jb=j-1;
      int ib=i-1;
      int ic=i+1;
      vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm) = 
           andiff(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),
                  v(k,j+offy_v-1,i+offx_v-1,icrm),irho(k,icrm))-
           (across(f(ind_f,k,jb+offy_s-1,ic+offx_s-1,icrm)+f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm)-
                   f(ind_f,k,jb+offy_s-1,ib+offx_s-1,icrm)-
                   f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),v(k,j+offy_v-1,i+offx_v-1,icrm), 
                   u(k,jb+offy_u-1,i+offx_u-1,icrm)+
                   u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+
                   u(k,jb+offy_u-1,ic+offx_u-1,icrm))+
            across(dd*(f(ind_f,kc,jb+offy_s-1,i+offx_s-1,icrm)+f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm)-
                   f(ind_f,kb,jb+offy_s-1,i+offx_s-1,icrm)-
                   f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm)),v(k,j+offy_v-1,i+offx_v-1,icrm), 
                   w(k,jb+offy_w-1,i+offx_w-1,icrm)+
                   w(k,j+offy_w-1,i+offx_w-1,icrm)+w(kc,j+offy_w-1,i+offx_w-1,icrm)+
                   w(kc,jb+offy_w-1,i+offx_w-1,icrm))) *irho(k,icrm);
    }
    if (i <= nx+1 && j <= ny+1) {
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      www(k,j+offy_www-1,i+offx_www-1,icrm) = 
           andiff(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),
                  w(k,j+offy_w-1,i+offx_w-1,icrm),irhow(k,icrm))-
          (across(f(ind_f,kb,j+offy_s-1,ic+offx_s-1,icrm)+f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm)-
                  f(ind_f,kb,j+offy_s-1,ib+offx_s-1,icrm)-
                  f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), 
                  u(kb,j+offy_u-1,i+offx_u-1,icrm)+
                  u(k,j+offy_u-1,i+offx_u-1,icrm)+u(k,j+offy_u-1,ic+offx_u-1,icrm)+
                  u(kb,j+offy_u-1,ic+offx_u-1,icrm))+
           across(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm)+f(ind_f,kb,jc+offy_s-1,i+offx_s-1,icrm)-
                  f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm)-
                  f(ind_f,kb,jb+offy_s-1,i+offx_s-1,icrm),w(k,j+offy_w-1,i+offx_w-1,icrm), 
                  v(kb,j+offy_v-1,i+offx_v-1,icrm)+
                  v(kb,jc+offy_v-1,i+offx_v-1,icrm)+v(k,jc+offy_v-1,i+offx_v-1,icrm)+
                  v(k,j+offy_v-1,i+offx_v-1,icrm))) *irho(k,icrm);
    }
  });

  //   for (int j=0; j<ny+4; j++) {
  //     for (int i=0; i<nx+4; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+4,nx+4,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int jb=j-1;
      int jc=j+1;
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm) = 
          max(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),max(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
          max(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),
          max(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),max(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),
          max(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
          max(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),mx(k,j,i,icrm))))))));
      mn(k,j,i,icrm) = 
          min(f(ind_f,k,j+offy_s-1,ib+offx_s-1,icrm),min(f(ind_f,k,j+offy_s-1,ic+offx_s-1,icrm),
          min(f(ind_f,k,jb+offy_s-1,i+offx_s-1,icrm),
          min(f(ind_f,k,jc+offy_s-1,i+offx_s-1,icrm),min(f(ind_f,kb,j+offy_s-1,i+offx_s-1,icrm),
          min(f(ind_f,kc,j+offy_s-1,i+offx_s-1,icrm),
          min(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm),mn(k,j,i,icrm))))))));
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+2; j++) {
    //     for (int i=0; i<nx+2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+2,nx+2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int jc=j+1;
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm))/
                ( pn3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pp3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pn3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pp3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pn3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pp3(www(k,j+offy_www-1,i+offx_www-1,icrm)))
                 *iadz(k,icrm)+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(ind_f,k,j+offy_s-1,i+offx_s-1,icrm)-mn(k,j,i,icrm))/
                ( pp3(uuu(k,j+offy_uuu-1,ic+offx_uuu-1,icrm)) + pn3(uuu(k,j+offy_uuu-1,i+offx_uuu-1,icrm))+
                  pp3(vvv(k,jc+offy_vvv-1,i+offx_vvv-1,icrm)) + pn3(vvv(k,j+offy_vvv-1,i+offx_vvv-1,icrm))+
                 (pp3(www(kc,j+offy_www-1,i+offx_www-1,icrm)) + pn3(www(k,j+offy_www-1,i+offx_www-1,icrm)))
                 *iadz(k,icrm)+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j <= ny-1) {
        int ib=i-1;
        uuu(k,j+offy_uuu,i+offx_uuu,icrm) = 
              pp3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), 
              mn(k,j+offy_m,ib+offx_m,icrm)))
             -pn3(uuu(k,j+offy_uuu,i+offx_uuu,icrm))*min(1.0,min(mx(k,j+offy_m,ib+offx_m,icrm),
             mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1) {
        int jb=j-1;
        vvv(k,j+offy_vvv,i+offx_vvv,icrm) =
              pp3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), 
              mn(k,jb+offy_m,i+offx_m,icrm)))
             -pn3(vvv(k,j+offy_vvv,i+offx_vvv,icrm))*min(1.0,min(mx(k,jb+offy_m,i+offx_m,icrm),
             mn(k,j+offy_m,i+offx_m,icrm)));
      }
      if (i <= nx-1 && j <= ny-1) {
        int kb=max(0,k-1);
        www(k,j+offy_www,i+offx_www,icrm) =
              pp3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(k,j+offy_m,i+offx_m,icrm), 
              mn(kb,j+offy_m,i+offx_m,icrm)))
             -pn3(www(k,j+offy_www,i+offx_www,icrm))*min(1.0,min(mx(kb,j+offy_m,i+offx_m,icrm),
             mn(k,j+offy_m,i+offx_m,icrm)));
        yakl::atomicAdd(flux(ind_flux,k,icrm),www(k,j+offy_www,i+offx_www,icrm));
      }
    });
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    int kc=k+1;
    f(ind_f,k,j+offy_s,i+offx_s,icrm) = 
         max(0.0,f(ind_f,k,j+offy_s,i+offx_s,icrm) -(uuu(k,j+offy_uuu,i+offx_uuu+1,icrm)-
                 uuu(k,j+offy_uuu,i+offx_uuu,icrm)+
                 vvv(k,j+offy_vvv+1,i+offx_vvv,icrm)-vvv(k,j+offy_vvv,i+offx_vvv,icrm)+
                 (www(k+1,j+offy_www,i+offx_www,icrm)-
                 www(k,j+offy_www,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}
