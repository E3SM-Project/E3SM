
#include "diffuse_mom2D.h"

void diffuse_mom2D(real5d &tk) {
  YAKL_SCOPE( dx            , :: dx );
  YAKL_SCOPE( dy            , :: dy );
  YAKL_SCOPE( dz            , :: dz );
  YAKL_SCOPE( adzw          , :: adzw );
  YAKL_SCOPE( grdf_x        , :: grdf_x );
  YAKL_SCOPE( grdf_z        , :: grdf_z );
  YAKL_SCOPE( u             , :: u );
  YAKL_SCOPE( v             , :: v );
  YAKL_SCOPE( w             , :: w );
  YAKL_SCOPE( na            , :: na );
  YAKL_SCOPE( dudt          , :: dudt );
  YAKL_SCOPE( dvdt          , :: dvdt );
  YAKL_SCOPE( dwdt          , :: dwdt );
  YAKL_SCOPE( uwsb          , :: uwsb );
  YAKL_SCOPE( vwsb          , :: vwsb );
  YAKL_SCOPE( rho           , :: rho );
  YAKL_SCOPE( rhow          , :: rhow );
  YAKL_SCOPE( fluxbu        , :: fluxbu );
  YAKL_SCOPE( fluxbv        , :: fluxbv );
  YAKL_SCOPE( fluxtu        , :: fluxtu );
  YAKL_SCOPE( fluxtv        , :: fluxtv );
  YAKL_SCOPE( adz           , :: adz );
  YAKL_SCOPE( ncrms         , :: ncrms );
  
  real4d fu("fu",nz,1,nx+1,ncrms);
  real4d fv("fv",nz,1,nx+1,ncrms);
  real4d fw("fw",nz,1,nx+1,ncrms);

  real rdx2=1.0/dx/dx;
  real rdx25=0.25*rdx2;
  int  constexpr j = 0;

  if (!docolumn) {
    // for (int k=0; k<nzm; k++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=k+1;
      int kcu=min(kc,nzm-1);
      real dxz=dx/(dz(icrm)*adzw(kc,icrm));
      real rdx21=rdx2    * grdf_x(k,icrm);
      real rdx251=rdx25  * grdf_x(k,icrm);
      int ic=i+1;
      real tkx=rdx21*tk(0,k,j+offy_d,i-1+offx_d,icrm);
      fu(k,j,i,icrm)=-2.0*tkx*(u(k,j+offy_u,ic-1+offx_u,icrm)-u(k,j+offy_u,i-1+offx_u,icrm));
      fv(k,j,i,icrm)=-tkx*(v(k,j+offy_v,ic+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm));
      tkx = rdx251*(tk(0,k,j+offy_d,i-1+offx_d,icrm)+tk(0,k,j+offy_d,ic-1+offx_d,icrm)+
            tk(0,kcu,j+offy_d,i-1+offx_d,icrm)+tk(0,kcu,j+offy_d,ic-1+offx_d,icrm));
      fw(k,j,i,icrm) = -tkx*(w(kc,j+offy_w,ic-1+offx_w,icrm)-w(kc,j+offy_w,i-1+offx_w,icrm)+
                       (u(kcu,j+offy_u,ic-1+offx_u,icrm)-u(k,j+offy_u,ic-1+offx_u,icrm))*dxz);
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kc=k+1;
      int ib=i-1;
      dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fu(k,j,i+1,icrm)-fu(k,j,ib+1,icrm));
      dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fv(k,j,i+1,icrm)-fv(k,j,ib+1,icrm));
      dwdt(na-1,kc,j,i,icrm)=dwdt(na-1,kc,j,i,icrm)-(fw(k,j,i+1,icrm)-fw(k,j,ib+1,icrm));
    });
  }

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    uwsb(k,icrm)=0.0;
    vwsb(k,icrm)=0.0;
  });

  // for (int k=0; k<nzm-1; k++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm-1,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=k+1;
    real rdz=1.0/dz(icrm);
    real rdz2 = rdz*rdz * grdf_z(k,icrm);
    real rdz25 = 0.25*rdz2;
    real dzx=dz(icrm)/dx;
    real iadz = 1.0/adz(k,icrm);
    real iadzw= 1.0/adzw(kc,icrm);
    int ib=i-1;
    real tkz=rdz2*tk(0,k,j+offy_d,i+offx_d,icrm);
    fw(kc,j,i+1,icrm)=-2.0*tkz*(w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rho(k,icrm)*iadz;
    tkz=rdz25*(tk(0,k,j+offy_d,i+offx_d,icrm)+tk(0,k,j+offy_d,ib+offx_d,icrm)+tk(0,kc,j+
        offy_d,i+offx_d,icrm)+tk(0,kc,j+offy_d,ib+offx_d,icrm));

    fu(kc,j,i+1,icrm)=-tkz*( (u(kc,j+offy_u,i+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*iadzw + 
                      (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*dzx)*rhow(kc,icrm);
    fv(kc,j,i+1,icrm)=-tkz*(v(kc,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*iadzw*rhow(kc,icrm);
    yakl::atomicAdd(uwsb(kc,icrm),fu(kc,j,i+1,icrm));
    yakl::atomicAdd(vwsb(kc,icrm),fv(kc,j,i+1,icrm));
   });
  
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nx,ncrms) , YAKL_LAMBDA (int i, int icrm) {
    real rdz=1.0/dz(icrm);
    real rdz2 = rdz*rdz * grdf_z(nzm-2,icrm);
    real tkz=rdz2*grdf_z(nzm-1,icrm)*tk(0,nzm-1,j+offy_d,i+offx_d,icrm);
    fw(nz-1,j,i+1,icrm)=-2.0*tkz*(w(nz-1,j+offy_d,i+offx_d,icrm)-w(nzm-1,j+offy_w,i+offx_w,icrm))/
                        adz(nzm-1,icrm)*rho(nzm-1,icrm);
    fu(0,j,i+1,icrm)=fluxbu(j,i,icrm) * rdz * rhow(0,icrm);
    fv(0,j,i+1,icrm)=fluxbv(j,i,icrm) * rdz * rhow(0,icrm);
    fu(nz-1,j,i+1,icrm)=fluxtu(j,i,icrm) * rdz * rhow(nz-1,icrm);
    fv(nz-1,j,i+1,icrm)=fluxtv(j,i,icrm) * rdz * rhow(nz-1,icrm);
    yakl::atomicAdd(uwsb(0,icrm),fu(0,j,i+1,icrm));
    yakl::atomicAdd(vwsb(0,icrm),fv(0,j,i+1,icrm));
  });
 
  // for (int k=0; k<nzm; k++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    int kc=k+1;
    real rhoi = 1.0/(rho(k,icrm)*adz(k,icrm));
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fu(kc,j,i+1,icrm)-fu(k,j,i+1,icrm))*rhoi;
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fv(kc,j,i+1,icrm)-fv(k,j,i+1,icrm))*rhoi;
  });

  // for (int k=0; k<nzm-1; k++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm-1,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
    real rhoi = 1.0/(rhow(k+1,icrm)*adzw(k+1,icrm));
    dwdt(na-1,k+1,j,i,icrm)=dwdt(na-1,k+1,j,i,icrm)-(fw(k+2,j,i+1,icrm)-fw(k+1,j,i+1,icrm))*rhoi;
  });

}

