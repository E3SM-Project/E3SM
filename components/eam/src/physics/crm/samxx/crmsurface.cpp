#include "crmsurface.h"

void crmsurface(real1d &bflx) {
  auto uhl     = ::uhl;
  auto vhl     = ::vhl;
  auto dtn     = ::dtn;
  auto utend   = ::utend;
  auto vtend   = ::vtend;
  auto taux0   = ::taux0;
  auto tauy0   = ::tauy0;
  auto u       = ::u;
  auto v       = ::v;
  auto rho     = ::rho;
  auto z       = ::z;
  auto z0      = ::z0;
  auto fluxbu  = ::fluxbu;
  auto fluxbv  = ::fluxbv;
  auto ug      = ::ug;
  auto vg      = ::vg;
  auto ncrms   = ::ncrms;

  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    uhl(icrm) = uhl(icrm) + dtn*utend(0,icrm);
    vhl(icrm) = vhl(icrm) + dtn*vtend(0,icrm);
    taux0(icrm) = 0.0;
    tauy0(icrm) = 0.0;
  });

  //  for (int j=0; j<ny; j++) {
  //   for (int i=0; i<nx; i++) {
  //     for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    real tmp2 = (0.5*(u(0,j+offy_u,i+1+offx_u,icrm)+u(0,j+offy_u,i+offx_u,icrm))+ug);
    real tmp3 = (0.5*(v(0,j+YES3D+offy_v,i+offx_v,icrm)+v(0,j+offy_v,i+offx_v,icrm))+vg);
    real u_h0 = max(1.0,sqrt(tmp2*tmp2+tmp3*tmp3));

    tmp2 = diag_ustar(z(0,icrm),bflx(icrm),u_h0,z0(icrm));
    real tau00 = rho(0,icrm)*tmp2*tmp2;

    fluxbu(j,i,icrm) = -(0.5*(u(0,j+offy_u,i+1+offx_u,icrm)+u(0,j+offy_u,i+offx_u,icrm))+ug-uhl(icrm))/u_h0*tau00;

    fluxbv(j,i,icrm) = -(0.5*(v(0,j+YES3D+offy_v,i+offx_v,icrm)+v(0,j+offy_v,i+offx_v,icrm))+vg-vhl(icrm))/u_h0*tau00;

    yakl::atomicAdd( taux0(icrm) , fluxbu(j,i,icrm)/( (real) nx * (real) ny ) );
    yakl::atomicAdd( tauy0(icrm) , fluxbv(j,i,icrm)/( (real) nx * (real) ny ) );
  });
}


