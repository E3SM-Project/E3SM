#include "advect_scalar2D.h"

void advect_scalar2D(real4d &f, real2d &flux) {
  auto &dowallx       = :: dowallx;
  auto &rank          = :: rank;
  auto &u             = :: u;
  auto &w             = :: w;
  auto &rho           = :: rho;
  auto &adz           = :: adz;
  auto &rhow          = :: rhow;
  auto &ncrms         = :: ncrms;

  bool constexpr nonos    = true;
  real constexpr eps      = 1.0e-10;
  int  constexpr offx_m   = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offx_www = 2;
  int  constexpr j        = 0;

  real4d mx   ("mx"   ,nzm,1,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,1,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,1,nx+5,ncrms);
  real4d www  ("www"  ,nz,1,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int i=0; i<nx+4; i++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<1-dimx1_u+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x==nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int iInd = i+ (nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(k,j,ib+offx_s-1,icrm),max(f(k,j,ic+offx_s-1,icrm),max(f(kb,j,i+offx_s-1,icrm),
                     max(f(kc,j,i+offx_s-1,icrm),f(k,j,i+offx_s-1,icrm)))));
      mn(k,j,i,icrm)=min(f(k,j,ib+offx_s-1,icrm),min(f(k,j,ic+offx_s-1,icrm),min(f(kb,j,i+offx_s-1,icrm),
                     min(f(kc,j,i+offx_s-1,icrm),f(k,j,i+offx_s-1,icrm)))));
    });
  }// nonos

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+5; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+5,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kb=max(0,k-1);
    uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(k,j,i-1+offx_s-2,icrm)+
                    min(0.0,u(k,j,i,icrm))*f(k,j,i+offx_s-2,icrm);
    if (i <= nx+3) {
      www(k,j,i,icrm)=max(0.0,w(k,j,i,icrm))*f(kb,j,i+offx_s-2,icrm)+min(0.0,w(k,j,i,icrm))*f(k,j,i+offx_s-2,icrm);
    }
    if (i == 1) {
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
  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+4,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    if (i >= 2 && i <= nx+1) {
      yakl::atomicAdd(flux(k,icrm),www(k,j,i,icrm));
    }
    f(k,j,i+offx_s-2,icrm) = f(k,j,i+offx_s-2,icrm) - (uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                                      (www(k+1,j,i,icrm)-www(k,j,i,icrm))*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+3; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+3,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=min(nzm-1,k+1);
    int kb=max(0,k-1);
    real dd=2.0/(kc-kb)/adz(k,icrm);
    int ib=i-1;
    uuu(k,j,i+offx_uuu-1,icrm) = 
         andiff2(f(k,j,ib+offx_s-1,icrm),f(k,j,i+offx_s-1,icrm),u(k,j,i+offx_u-1,icrm),irho(k,icrm)) - 
         across2(dd*(f(kc,j,ib+offx_s-1,icrm)+f(kc,j,i+offx_s-1,icrm)-f(kb,j,ib+offx_s-1,icrm)-f(kb,j,i+offx_s-1,icrm)),
                 u(k,j,i+offx_u-1,icrm), w(k,j,ib+offx_w-1,icrm)+w(kc,j,ib+offx_w-1,icrm)+w(k,j,i+offx_w-1,icrm)+
                 w(kc,j,i+offx_w-1,icrm)) *irho(k,icrm);
    if (i <= nxp1) {
      int ic=i+1;
      www(k,j,i+offx_www-1,icrm) = 
         andiff2(f(kb,j,i+offx_s-1,icrm),f(k,j,i+offx_s-1,icrm),w(k,j,i+offx_w-1,icrm),irhow(k,icrm)) -
         across2(f(kb,j,ic+offx_s-1,icrm)+f(k,j,ic+offx_s-1,icrm)-f(kb,j,ib+offx_s-1,icrm)-f(k,j,ib+offx_s-1,icrm),
                 w(k,j,i+offx_w-1,icrm), u(kb,j,i+offx_u-1,icrm)+u(k,j,i+offx_u-1,icrm)+u(k,j,ic+offx_u-1,icrm)+
                 u(kb,j,ic+offx_u-1,icrm)) *irho(k,icrm);
    }
  });

  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(k,j,ib+offx_s-1,icrm),max(f(k,j,ic+offx_s-1,icrm),max(f(kb,j,i+offx_s-1,icrm),
                     max(f(kc,j,i+offx_s-1,icrm),max(f(k,j,i+offx_s-1,icrm),mx(k,j,i,icrm))))));
      mn(k,j,i,icrm)=min(f(k,j,ib+offx_s-1,icrm),min(f(k,j,ic+offx_s-1,icrm),min(f(kb,j,i+offx_s-1,icrm),
                     min(f(kc,j,i+offx_s-1,icrm),min(f(k,j,i+offx_s-1,icrm),mn(k,j,i,icrm))))));
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(k,j,i+offx_s-1,icrm))/(pn2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pp2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pn2(www(kc,j,i+offx_www-1,icrm)) +
                     pp2(www(k,j,i+offx_www-1,icrm)))+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(k,j,i+offx_s-1,icrm)-mn(k,j,i,icrm))/(pp2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pn2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pp2(www(kc,j,i+offx_www-1,icrm)) +
                     pn2(www(k,j,i+offx_www-1,icrm)))+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+1; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int ib=i-1;
      uuu(k,j,i+offx_uuu,icrm) =
            pp2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(k,j,ib+offx_m,icrm))) -
            pn2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,ib+offx_m,icrm),mn(k,j,i+offx_m,icrm)));
      if (i <= nx-1) {
        int kb=max(0,k-1);
        www(k,j,i+offx_www,icrm) =
            pp2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(kb,j,i+offx_m,icrm))) -
            pn2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(kb,j,i+offx_m,icrm),mn(k,j,i+offx_m,icrm)));
        yakl::atomicAdd(flux(k,icrm), www(k,j,i+offx_www,icrm));
      }
    });
  } // nonos

  // for (int k=0; k<nzm; k++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=k+1;
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    f(k,j,i+offx_s,icrm)= max(0.0, f(k,j,i+offx_s,icrm) - (uuu(k,j,i+1+offx_uuu,icrm)-uuu(k,j,i+offx_uuu,icrm) +
                         (www(k+1,j,i+offx_www,icrm)-www(k,j,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}



void advect_scalar2D(real5d &f, int ind_f, real2d &flux) {
  auto &dowallx       = :: dowallx;
  auto &rank          = :: rank;
  auto &u             = :: u;
  auto &w             = :: w;
  auto &rho           = :: rho;
  auto &adz           = :: adz;
  auto &rhow          = :: rhow;
  auto &ncrms         = :: ncrms;

  bool constexpr nonos = true;
  real constexpr eps = 1.0e-10;
  int  constexpr offx_m = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offx_www = 2;
  int  constexpr j = 0;

  real4d mx   ("mx"   ,nzm,1,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,1,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,1,nx+5,ncrms);
  real4d www  ("www"  ,nz,1,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int i=0; i<nx+4; i++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<1-dimx1_u+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x==nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int iInd = i+ (nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(ind_f,k,j,ib+offx_s-1,icrm),max(f(ind_f,k,j,ic+offx_s-1,icrm),
                     max(f(ind_f,kb,j,i+offx_s-1,icrm),max(f(ind_f,kc,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm)))));
      mn(k,j,i,icrm)=min(f(ind_f,k,j,ib+offx_s-1,icrm),min(f(ind_f,k,j,ic+offx_s-1,icrm),
                     min(f(ind_f,kb,j,i+offx_s-1,icrm),min(f(ind_f,kc,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm)))));
    });
  }// nonos

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+5; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+5,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kb=max(0,k-1);
    uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(ind_f,k,j,i-1+offx_s-2,icrm)+
                    min(0.0,u(k,j,i,icrm))*f(ind_f,k,j,i+offx_s-2,icrm);
    if (i <= nx+3) {
      www(k,j,i,icrm) = max(0.0,w(k,j,i,icrm))*f(ind_f,kb,j,i+offx_s-2,icrm)+
                        min(0.0,w(k,j,i,icrm))*f(ind_f,k,j,i+offx_s-2,icrm);
    }
    if (i == 1) {
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
  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+4,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    if (i >= 2 && i <= nx+1) {
      yakl::atomicAdd(flux(k,icrm),www(k,j,i,icrm));
    }
    f(ind_f,k,j,i+offx_s-2,icrm) = f(ind_f,k,j,i+offx_s-2,icrm) -
                                   (uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                   (www(k+1,j,i,icrm)-www(k,j,i,icrm))*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+3; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+3,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=min(nzm-1,k+1);
    int kb=max(0,k-1);
    real dd=2.0/(kc-kb)/adz(k,icrm);
    int ib=i-1;
    uuu(k,j,i+offx_uuu-1,icrm) = 
          andiff2(f(ind_f,k,j,ib+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm),u(k,j,i+offx_u-1,icrm),irho(k,icrm)) - 
          across2(dd*(f(ind_f,kc,j,ib+offx_s-1,icrm)+f(ind_f,kc,j,i+offx_s-1,icrm)-
                  f(ind_f,kb,j,ib+offx_s-1,icrm)-f(ind_f,kb,j,i+offx_s-1,icrm)),
                  u(k,j,i+offx_u-1,icrm), w(k,j,ib+offx_w-1,icrm)+w(kc,j,ib+offx_w-1,icrm)+
                  w(k,j,i+offx_w-1,icrm)+w(kc,j,i+offx_w-1,icrm)) *irho(k,icrm);
    if (i <= nxp1) {
      int ic=i+1;
      www(k,j,i+offx_www-1,icrm) = 
          andiff2(f(ind_f,kb,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm),w(k,j,i+offx_w-1,icrm),irhow(k,icrm)) - 
          across2(f(ind_f,kb,j,ic+offx_s-1,icrm)+f(ind_f,k,j,ic+offx_s-1,icrm)-
                  f(ind_f,kb,j,ib+offx_s-1,icrm)-f(ind_f,k,j,ib+offx_s-1,icrm),
                  w(k,j,i+offx_w-1,icrm), u(kb,j,i+offx_u-1,icrm)+u(k,j,i+offx_u-1,icrm)+
                  u(k,j,ic+offx_u-1,icrm)+u(kb,j,ic+offx_u-1,icrm)) *irho(k,icrm);
    }
  });

  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(ind_f,k,j,ib+offx_s-1,icrm),max(f(ind_f,k,j,ic+offx_s-1,icrm),max(f(ind_f,kb,j,i+offx_s-1,icrm),
                     max(f(ind_f,kc,j,i+offx_s-1,icrm),max(f(ind_f,k,j,i+offx_s-1,icrm),mx(k,j,i,icrm))))));
      mn(k,j,i,icrm)=min(f(ind_f,k,j,ib+offx_s-1,icrm),min(f(ind_f,k,j,ic+offx_s-1,icrm),min(f(ind_f,kb,j,i+offx_s-1,icrm),
                     min(f(ind_f,kc,j,i+offx_s-1,icrm),min(f(ind_f,k,j,i+offx_s-1,icrm),mn(k,j,i,icrm))))));
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(ind_f,k,j,i+offx_s-1,icrm))/(pn2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pp2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pn2(www(kc,j,i+offx_www-1,icrm)) +
                     pp2(www(k,j,i+offx_www-1,icrm)))+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(ind_f,k,j,i+offx_s-1,icrm)-mn(k,j,i,icrm))/(pp2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pn2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pp2(www(kc,j,i+offx_www-1,icrm)) +
                     pn2(www(k,j,i+offx_www-1,icrm)))+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+1; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int ib=i-1;
      uuu(k,j,i+offx_uuu,icrm)= pp2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(k,j,ib+offx_m,icrm))) -
                       pn2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,ib+offx_m,icrm),mn(k,j,i+offx_m,icrm)));
      if (i <= nx-1) {
        int kb=max(0,k-1);
        www(k,j,i+offx_www,icrm)= pp2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(kb,j,i+offx_m,icrm))) -
                         pn2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(kb,j,i+offx_m,icrm),mn(k,j,i+offx_m,icrm)));

        yakl::atomicAdd(flux(k,icrm), www(k,j,i+offx_www,icrm));
      }
    });
  } // nonos

  // for (int k=0; k<nzm; k++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=k+1;
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    f(ind_f,k,j,i+offx_s,icrm)= max(0.0, f(ind_f,k,j,i+offx_s,icrm) - (uuu(k,j,i+1+offx_uuu,icrm)-uuu(k,j,i+offx_uuu,icrm) +
                   (www(k+1,j,i+offx_www,icrm)-www(k,j,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}

void advect_scalar2D(real5d &f, int ind_f, real3d &flux, int ind_flux) {
  auto &dowallx       = :: dowallx;
  auto &rank          = :: rank;
  auto &u             = :: u;
  auto &w             = :: w;
  auto &rho           = :: rho;
  auto &adz           = :: adz;
  auto &rhow          = :: rhow;
  auto &ncrms         = :: ncrms;

  bool constexpr nonos = true;
  real constexpr eps = 1.0e-10;
  int  constexpr offx_m = 1;
  int  constexpr offx_uuu = 2;
  int  constexpr offx_www = 2;
  int  constexpr j = 0;

  real4d mx   ("mx"   ,nzm,1,nx+2,ncrms);
  real4d mn   ("mn"   ,nzm,1,nx+2,ncrms);
  real4d uuu  ("uuu"  ,nzm,1,nx+5,ncrms);
  real4d www  ("www"  ,nz,1,nx+4,ncrms);
  real2d iadz ("iadz" ,nzm,ncrms);
  real2d irho ("irho" ,nzm,ncrms);
  real2d irhow("irhow",nzm,ncrms);

  // for (int i=0; i<nx+4; i++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(nz-1,j,i,icrm)=0.0;
  });

  if (dowallx) {
    if (rank%nsubdomains_x == 0) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<1-dimx1_u+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        u(k,j,i,icrm) = 0.0;
      });
    }
    if (rank%nsubdomains_x==nsubdomains_x-1) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<dimx2_u-(nx+1)+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int iInd = i+ (nx+2);
        u(k,j,iInd,icrm) = 0.0;
      });
    }
  }

  if (nonos) {
    
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(ind_f,k,j,ib+offx_s-1,icrm),max(f(ind_f,k,j,ic+offx_s-1,icrm),max(f(ind_f,kb,j,i+offx_s-1,icrm),
                     max(f(ind_f,kc,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm)))));
      mn(k,j,i,icrm)=min(f(ind_f,k,j,ib+offx_s-1,icrm),min(f(ind_f,k,j,ic+offx_s-1,icrm),min(f(ind_f,kb,j,i+offx_s-1,icrm),
                     min(f(ind_f,kc,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm)))));
    });
  }// nonos

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+5; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+5,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kb=max(0,k-1);
    uuu(k,j,i,icrm)=max(0.0,u(k,j,i,icrm))*f(ind_f,k,j,i-1+offx_s-2,icrm)+
                    min(0.0,u(k,j,i,icrm))*f(ind_f,k,j,i+offx_s-2,icrm);
    if (i <= nx+3) {
      www(k,j,i,icrm)=max(0.0,w(k,j,i,icrm))*f(ind_f,kb,j,i+offx_s-2,icrm)+min(0.0,w(k,j,i,icrm))*f(ind_f,k,j,i+offx_s-2,icrm);
    }
    if (i == 1) {
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
  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+4,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    if (i >= 2 && i <= nx+1) {
      yakl::atomicAdd(flux(ind_flux,k,icrm),www(k,j,i,icrm));
    }
    f(ind_f,k,j,i+offx_s-2,icrm) = f(ind_f,k,j,i+offx_s-2,icrm) - (uuu(k,j,i+1,icrm)-uuu(k,j,i,icrm) +
                                   (www(k+1,j,i,icrm)-www(k,j,i,icrm))*iadz(k,icrm))*irho(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx+3; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx+3,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=min(nzm-1,k+1);
    int kb=max(0,k-1);
    real dd=2.0/(kc-kb)/adz(k,icrm);
    int ib=i-1;
    uuu(k,j,i+offx_uuu-1,icrm) = 
         andiff2(f(ind_f,k,j,ib+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm),u(k,j,i+offx_u-1,icrm),irho(k,icrm)) - 
         across2(dd*(f(ind_f,kc,j,ib+offx_s-1,icrm)+f(ind_f,kc,j,i+offx_s-1,icrm)-
                 f(ind_f,kb,j,ib+offx_s-1,icrm)-f(ind_f,kb,j,i+offx_s-1,icrm)),
                 u(k,j,i+offx_u-1,icrm), w(k,j,ib+offx_w-1,icrm)+w(kc,j,ib+offx_w-1,icrm)+
                 w(k,j,i+offx_w-1,icrm)+w(kc,j,i+offx_w-1,icrm)) *irho(k,icrm);
    if (i <= nxp1) {
      int ic=i+1;
      www(k,j,i+offx_www-1,icrm) = 
         andiff2(f(ind_f,kb,j,i+offx_s-1,icrm),f(ind_f,k,j,i+offx_s-1,icrm),w(k,j,i+offx_w-1,icrm),irhow(k,icrm)) - 
         across2(f(ind_f,kb,j,ic+offx_s-1,icrm)+f(ind_f,k,j,ic+offx_s-1,icrm)-
                 f(ind_f,kb,j,ib+offx_s-1,icrm)-f(ind_f,k,j,ib+offx_s-1,icrm),
                 w(k,j,i+offx_w-1,icrm), u(kb,j,i+offx_u-1,icrm)+u(k,j,i+offx_u-1,icrm)+
                 u(k,j,ic+offx_u-1,icrm)+u(kb,j,ic+offx_u-1,icrm)) *irho(k,icrm);
    }
  });

  //  for (int i=0; i<nx+4; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx+4,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    www(0,j,i,icrm) = 0.0;
  });

  if (nonos) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int kb=max(0,k-1);
      int ib=i-1;
      int ic=i+1;
      mx(k,j,i,icrm)=max(f(ind_f,k,j,ib+offx_s-1,icrm),max(f(ind_f,k,j,ic+offx_s-1,icrm),max(f(ind_f,kb,j,i+offx_s-1,icrm),
                     max(f(ind_f,kc,j,i+offx_s-1,icrm),max(f(ind_f,k,j,i+offx_s-1,icrm),mx(k,j,i,icrm))))));
      mn(k,j,i,icrm)=min(f(ind_f,k,j,ib+offx_s-1,icrm),min(f(ind_f,k,j,ic+offx_s-1,icrm),min(f(ind_f,kb,j,i+offx_s-1,icrm),
                     min(f(ind_f,kc,j,i+offx_s-1,icrm),min(f(ind_f,k,j,i+offx_s-1,icrm),mn(k,j,i,icrm))))));
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+2; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+2,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=min(nzm-1,k+1);
      int ic=i+1;
      mx(k,j,i,icrm)=rho(k,icrm)*(mx(k,j,i,icrm)-f(ind_f,k,j,i+offx_s-1,icrm))/(pn2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pp2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pn2(www(kc,j,i+offx_www-1,icrm)) +
                     pp2(www(k,j,i+offx_www-1,icrm)))+eps);
      mn(k,j,i,icrm)=rho(k,icrm)*(f(ind_f,k,j,i+offx_s-1,icrm)-mn(k,j,i,icrm))/(pp2(uuu(k,j,ic+offx_uuu-1,icrm)) +
                     pn2(uuu(k,j,i+offx_uuu-1,icrm))+iadz(k,icrm)*(pp2(www(kc,j,i+offx_www-1,icrm)) +
                     pn2(www(k,j,i+offx_www-1,icrm)))+eps);
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx+1; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int ib=i-1;
      uuu(k,j,i+offx_uuu,icrm)= pp2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(k,j,ib+offx_m,icrm))) -
                                pn2(uuu(k,j,i+offx_uuu,icrm))*min(1.0,min(mx(k,j,ib+offx_m,icrm),mn(k,j,i+offx_m,icrm)));
      if (i <= nx-1) {
        int kb=max(0,k-1);
        www(k,j,i+offx_www,icrm)= pp2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(k,j,i+offx_m,icrm), mn(kb,j,i+offx_m,icrm))) -
                                  pn2(www(k,j,i+offx_www,icrm))*min(1.0,min(mx(kb,j,i+offx_m,icrm),mn(k,j,i+offx_m,icrm)));

        yakl::atomicAdd(flux(ind_flux,k,icrm), www(k,j,i+offx_www,icrm));
      }
    });
  } // nonos

  // for (int k=0; k<nzm; k++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=k+1;
    // MK: added fix for very small negative values (relative to positive values)
    //     especially  when such large numbers as
    //     hydrometeor concentrations are advected. The reason for negative values is
    //     most likely truncation error.
    f(ind_f,k,j,i+offx_s,icrm)= max(0.0, f(ind_f,k,j,i+offx_s,icrm) - (uuu(k,j,i+1+offx_uuu,icrm)-uuu(k,j,i+offx_uuu,icrm) +
                                (www(k+1,j,i+offx_www,icrm)-www(k,j,i+offx_www,icrm))*iadz(k,icrm))*irho(k,icrm));
  });

}
