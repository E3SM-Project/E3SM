
#include "damping.h"

void damping() {
  YAKL_SCOPE( z              , ::z );
  YAKL_SCOPE( u              , ::u );
  YAKL_SCOPE( v              , ::v );
  YAKL_SCOPE( t              , ::t );
  YAKL_SCOPE( na             , ::na );
  YAKL_SCOPE( dudt           , ::dudt );
  YAKL_SCOPE( dvdt           , ::dvdt );
  YAKL_SCOPE( dwdt           , ::dwdt );
  YAKL_SCOPE( w              , ::w );
  YAKL_SCOPE( dtn            , ::dtn );
  YAKL_SCOPE( micro_field    , ::micro_field );
  YAKL_SCOPE( qv             , ::qv );
  YAKL_SCOPE( qv0            , ::qv0 );
  YAKL_SCOPE( ncrms          , ::ncrms );

  real constexpr tau_min    = 60.0;
  real constexpr tau_max    = 450.0;
  real constexpr fractional_damp_depth = 0.4;

  int1d  n_damp    ("n_damp",ncrms);
  int2d  do_damping("n_damp",nzm,ncrms);
  real2d t0loc     ("t0loc" ,nzm,ncrms);
  real2d u0loc     ("u0loc" ,nzm,ncrms);
  real2d v0loc     ("v0loc" ,nzm,ncrms);
  real2d tau       ("tau"   ,nzm,ncrms);

  if (tau_min < 2.0*dt) { 
    std::cout << "Error: in damping() tau_min is too small!";
    exit(-1);
  }

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    n_damp(icrm) = 0;
  });

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    if(z(nzm-1,icrm)-z(k,icrm) < fractional_damp_depth*z(nzm-1,icrm)) {
      do_damping(k,icrm)=1;
    } else {
      do_damping(k,icrm)=0;
    }
    yakl::atomicAdd( n_damp(icrm), do_damping(k,icrm) );
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    tau(k,icrm) = 0;
    if ( (k <= nzm-1) && (k >= nzm-1-n_damp(icrm)) ) {
      tau(k,icrm) = tau_min * pow( (tau_max/tau_min) ,
                                 ( ( z(nzm-1,icrm) - z(k,icrm) ) / ( z(nzm-1,icrm) - z( nzm-1-n_damp(icrm) , icrm ) ) ) );
      tau(k,icrm) = 1. / tau(k,icrm);
    }
  });

  // recalculate grid-mean u0, v0, t0 first,
  // as t has been updated. No need for qv0, as
  // qv has not been updated yet the calculation of qv0.

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    u0loc(k,icrm)=0.0;
    v0loc(k,icrm)=0.0;
    t0loc(k,icrm)=0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real tmp;

    tmp = u(k,offy_u+j,offx_u+i,icrm)/( (real) nx * (real) ny );
    yakl::atomicAdd(u0loc(k,icrm),tmp);

    tmp = v(k,offy_v+j,offx_v+i,icrm)/( (real) nx * (real) ny );
    yakl::atomicAdd(v0loc(k,icrm),tmp);

    tmp = t(k,offy_s+j,offx_s+i,icrm)/( (real) nx * (real) ny );
    yakl::atomicAdd(t0loc(k,icrm),tmp);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int idwv = index_water_vapor;
    if ( k <= nzm-1 && k >= nzm-1-n_damp(icrm) ) {
      dudt       (na-1,k,       j,       i,icrm) -=     (u (k,offy_u+j,offx_u+i,icrm)-u0loc(k,icrm)) * tau(k,icrm);
      dvdt       (na-1,k,       j,       i,icrm) -=     (v (k,offy_v+j,offx_v+i,icrm)-v0loc(k,icrm)) * tau(k,icrm);
      dwdt       (na-1,k,       j,       i,icrm) -=      w (k,offy_w+j,offx_w+i,icrm)                * tau(k,icrm);
      t          (     k,offy_s+j,offx_s+i,icrm) -= dtn*(t (k,offy_s+j,offx_s+i,icrm)-t0loc(k,icrm)) * tau(k,icrm);
      micro_field(idwv,k,offy_s+j,offx_s+i,icrm) -= dtn*(qv(k,       j,       i,icrm)-qv0  (k,icrm)) * tau(k,icrm);
    }
  });

}

