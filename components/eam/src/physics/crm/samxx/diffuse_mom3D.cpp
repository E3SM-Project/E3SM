#include "diffuse_mom3D.h"

void diffuse_mom3D(real5d &tk) {
  auto &dx            = :: dx;
  auto &dy            = :: dy;
  auto &dz            = :: dz;
  auto &adzw          = :: adzw;
  auto &grdf_x        = :: grdf_x;
  auto &grdf_y        = :: grdf_y;
  auto &grdf_z        = :: grdf_z;
  auto &u             = :: u;
  auto &v             = :: v;
  auto &w             = :: w;
  auto &na            = :: na;
  auto &dudt          = :: dudt;
  auto &dvdt          = :: dvdt;
  auto &dwdt          = :: dwdt;
  auto &uwsb          = :: uwsb;
  auto &vwsb          = :: vwsb;
  auto &rho           = :: rho;
  auto &rhow          = :: rhow;
  auto &fluxbu        = :: fluxbu;
  auto &fluxbv        = :: fluxbv;
  auto &fluxtu        = :: fluxtu;
  auto &fluxtv        = :: fluxtv;
  auto &adz           = :: adz;
  auto &ncrms         = :: ncrms;

  real4d fu("fu",nz,ny+1,nx+1,ncrms);
  real4d fv("fv",nz,ny+1,nx+1,ncrms);
  real4d fw("fw",nz,ny+1,nx+1,ncrms);

  real rdx2=1.0/(dx*dx);
  real rdy2=1.0/(dy*dy);
  real rdx25=0.25*rdx2;
  real rdy25=0.25*rdy2;
  real dxy=dx/dy;
  real dyx=dy/dx;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx+1; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int jb=j-1;
    int kc=k+1;
    int kcu=min(kc,nzm-1);
    int ic=i+1;
    real dxz=dx/(dz(icrm)*adzw(kc,icrm));
    real rdx21=rdx2    * grdf_x(k,icrm);
    real rdx251=rdx25  * grdf_x(k,icrm);
    real tkx=rdx21*tk(0,k,j+offy_d,i-1+offx_d,icrm);
    fu(k,j+1,i,icrm)=-2.0*tkx*(u(k,j+offy_u,ic-1+offx_u,icrm)-u(k,j+offy_u,i-1+offx_u,icrm));
    tkx=rdx251*(tk(0,k,j+offy_d,i-1+offx_d,icrm)+tk(0,k,jb+offy_d,i-1+offx_d,icrm)+
                tk(0,k,j+offy_d,ic-1+offx_d,icrm)+tk(0,k,jb+offy_d,ic-1+offx_d,icrm));
    fv(k,j+1,i,icrm)=-tkx*(v(k,j+offy_v,ic-1+offx_v,icrm)-v(k,j+offy_v,i-1+offx_v,icrm)+
                          (u(k,j+offy_u,ic-1+offx_u,icrm)-u(k,jb+offy_u,ic-1+offx_u,icrm))*dxy);
    tkx=rdx251*(tk(0,k,j+offy_d,i-1+offx_d,icrm)+tk(0,k,j+offy_d,ic-1+offx_d,icrm)+
                tk(0,kcu,j+offy_d,i-1+offx_d,icrm)+tk(0,kcu,j+offy_d,ic-1+offx_d,icrm));
    fw(k,j+1,i,icrm)=-tkx*(w(kc,j+offy_w,ic-1+offx_w,icrm)-w(kc,j+offy_w,i-1+offx_w,icrm)+
                          (u(kcu,j+offy_u,ic-1+offx_u,icrm)-u(k,j+offy_u,ic-1+offx_u,icrm))*dxz);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kc=k+1;
    int ib=i-1;
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fu(k,j+1,i+1,icrm)-fu(k,j+1,ib+1,icrm));
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fv(k,j+1,i+1,icrm)-fv(k,j+1,ib+1,icrm));
    dwdt(na-1,kc,j,i,icrm)=dwdt(na-1,kc,j,i,icrm)-(fw(k,j+1,i+1,icrm)-fw(k,j+1,ib+1,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny+1; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny+1,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int jc=j+1;
    int kc=k+1;
    int kcu=min(kc,nzm-1);
    int ib=i-1;
    real dyz=dy/(dz(icrm)*adzw(kc,icrm));
    real rdy21=rdy2    * grdf_y(k,icrm);
    real rdy251=rdy25  * grdf_y(k,icrm);
    real tky=rdy21*tk(0,k,j-1+offy_d,i+offx_d,icrm);
    fv(k,j,i+1,icrm)=-2.0*tky*(v(k,jc-1+offy_v,i+offx_v,icrm)-v(k,j-1+offy_v,i+offx_v,icrm));
    tky=rdy251*(tk(0,k,j-1+offy_d,i+offx_d,icrm)+tk(0,k,j-1+offy_d,ib+offx_d,icrm)+
                tk(0,k,jc-1+offy_d,i+offx_d,icrm)+tk(0,k,jc-1+offy_d,ib+offx_d,icrm));
    fu(k,j,i+1,icrm)=-tky*(u(k,jc-1+offy_u,i+offx_u,icrm)-u(k,j-1+offy_u,i+offx_u,icrm)+
                          (v(k,jc-1+offy_v,i+offx_v,icrm)-v(k,jc-1+offy_v,ib+offx_v,icrm))*dyx);
    tky=rdy251*(tk(0,k,j-1+offy_d,i+offx_d,icrm)+tk(0,k,jc-1+offy_d,i+offx_d,icrm)+
                tk(0,kcu,j-1+offy_d,i+offx_d,icrm)+tk(0,kcu,jc-1+offy_d,i+offx_d,icrm));
    fw(k,j,i+1,icrm)=-tky*(w(kc,jc-1+offy_w,i+offx_w,icrm)-w(kc,j-1+offy_w,i+offx_w,icrm)+
                          (v(kcu,jc-1+offy_v,i+offx_v,icrm)-v(k,jc-1+offy_v,i+offx_v,icrm))*dyz);
  });


  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int jb=j-1;
    int kc=k+1;
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fu(k,j+1,i+1,icrm)-fu(k,jb+1,i+1,icrm));
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fv(k,j+1,i+1,icrm)-fv(k,jb+1,i+1,icrm));
    dwdt(na-1,kc,j,i,icrm)=dwdt(na-1,kc,j,i,icrm)-(fw(k,j+1,i+1,icrm)-fw(k,jb+1,i+1,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    uwsb(k,icrm)=0.0;
    vwsb(k,icrm)=0.0;
  });

  // for (int k=0; k<nzm-1; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int jb=j-1;
    int kc=k+1;
    int ib=i-1;
    real rdz=1.0/dz(icrm);
    real rdz2 = rdz*rdz * grdf_z(k,icrm);
    real rdz25 = 0.25*rdz2;
    real iadz = 1.0/adz(k,icrm);
    real iadzw= 1.0/adzw(kc,icrm);
    real dzx=dz(icrm)/dx;
    real dzy=dz(icrm)/dy;
    real tkz=rdz2*tk(0,k,j+offy_d,i+offx_d,icrm);
    fw(kc,j+1,i+1,icrm)=-2.0*tkz*(w(kc,j+offy_w,i+offx_w,icrm)-w(k,j+offy_w,i+offx_w,icrm))*rho(k,icrm)*iadz;
    tkz=rdz25*(tk(0,k,j+offy_d,i+offx_d,icrm)+tk(0,k,j+offy_d,ib+offx_d,icrm)+tk(0,kc,j+offy_d,i+offx_d,icrm)+
               tk(0,kc,j+offy_d,ib+offx_d,icrm));
    fu(kc,j+1,i+1,icrm)=-tkz*( (u(kc,j+offy_u,i+offx_u,icrm)-u(k,j+offy_u,i+offx_u,icrm))*iadzw + 
                               (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,j+offy_w,ib+offx_w,icrm))*dzx)*rhow(kc,icrm);
    tkz=rdz25*(tk(0,k,j+offy_d,i+offx_d,icrm)+tk(0,k,jb+offy_d,i+offx_d,icrm)+tk(0,kc,j+offy_d,i+offx_d,icrm)+
               tk(0,kc,jb+offy_d,i+offx_d,icrm));
    fv(kc,j+1,i+1,icrm)=-tkz*( (v(kc,j+offy_v,i+offx_v,icrm)-v(k,j+offy_v,i+offx_v,icrm))*iadzw + 
                               (w(kc,j+offy_w,i+offx_w,icrm)-w(kc,jb+offy_w,i+offx_w,icrm))*dzy)*rhow(kc,icrm);
    yakl::atomicAdd(uwsb(kc,icrm),fu(kc,j+1,i+1,icrm));
    yakl::atomicAdd(vwsb(kc,icrm),fv(kc,j+1,i+1,icrm));
  });

  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    real rdz=1.0/dz(icrm);
    real rdz2 = rdz*rdz * grdf_z(nzm-2,icrm);
    real tkz=rdz2*grdf_z(nzm-1,icrm)*tk(0,nzm-1,j+offy_d,i+offx_d,icrm);
    fw(nz-1,j+1,i+1,icrm)=-2.0*tkz*(w(nz-1,j+offy_d,i+offx_d,icrm)-w(nzm-1,j+offy_w,i+offx_w,icrm))/
                          adz(nzm-1,icrm)*rho(nzm-1,icrm);
    fu(0,j+1,i+1,icrm)=fluxbu(j,i,icrm) * rdz * rhow(0,icrm);
    fv(0,j+1,i+1,icrm)=fluxbv(j,i,icrm) * rdz * rhow(0,icrm);
    fu(nz-1,j+1,i+1,icrm)=fluxtu(j,i,icrm) * rdz * rhow(nz-1,icrm);
    fv(nz-1,j+1,i+1,icrm)=fluxtv(j,i,icrm) * rdz * rhow(nz-1,icrm);
    yakl::atomicAdd(uwsb(0,icrm),fu(0,j+1,i+1,icrm));
    yakl::atomicAdd(vwsb(0,icrm),fv(0,j+1,i+1,icrm));
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int kc=k+1;
    real rhoi = 1.0/(rho(k,icrm)*adz(k,icrm));
    dudt(na-1,k,j,i,icrm)=dudt(na-1,k,j,i,icrm)-(fu(kc,j+1,i+1,icrm)-fu(k,j+1,i+1,icrm))*rhoi;
    dvdt(na-1,k,j,i,icrm)=dvdt(na-1,k,j,i,icrm)-(fv(kc,j+1,i+1,icrm)-fv(k,j+1,i+1,icrm))*rhoi;
  });

  // for (int k=0; k<nzm-1; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real rhoi = 1.0/(rhow(k+1,icrm)*adzw(k+1,icrm));
    dwdt(na-1,k+1,j,i,icrm)=dwdt(na-1,k+1,j,i,icrm)-(fw(k+2,j+1,i+1,icrm)-fw(k+1,j+1,i+1,icrm))*rhoi;
  });

}
