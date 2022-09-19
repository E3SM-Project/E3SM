#include "crmsurface.h"

void crmsurface(real1d &bflx) {
  YAKL_SCOPE( uhl      , ::uhl);
  YAKL_SCOPE( vhl      , ::vhl);
  YAKL_SCOPE( dtn      , ::dtn);
  YAKL_SCOPE( utend    , ::utend);
  YAKL_SCOPE( vtend    , ::vtend);
  YAKL_SCOPE( taux0    , ::taux0);
  YAKL_SCOPE( tauy0    , ::tauy0);
  YAKL_SCOPE( u        , ::u);
  YAKL_SCOPE( v        , ::v);
  YAKL_SCOPE( rho      , ::rho);
  YAKL_SCOPE( z        , ::z);
  YAKL_SCOPE( z0       , ::z0);
  YAKL_SCOPE( fluxbu   , ::fluxbu);
  YAKL_SCOPE( fluxbv   , ::fluxbv);
  YAKL_SCOPE( ug       , ::ug);
  YAKL_SCOPE( vg       , ::vg);
  YAKL_SCOPE( ncrms    , ::ncrms);

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


