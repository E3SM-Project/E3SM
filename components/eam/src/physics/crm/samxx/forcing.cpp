
#include "forcing.h"

void forcing() {
  YAKL_SCOPE( ncrms         , ::ncrms );
  YAKL_SCOPE( t             , ::t );
  YAKL_SCOPE( ttend         , ::ttend );
  YAKL_SCOPE( dtn           , ::dtn );
  YAKL_SCOPE( micro_field   , ::micro_field );
  YAKL_SCOPE( qtend         , ::qtend );
  YAKL_SCOPE( dudt          , ::dudt );
  YAKL_SCOPE( dvdt          , ::dvdt );
  YAKL_SCOPE( na            , ::na );
  YAKL_SCOPE( utend         , ::utend );
  YAKL_SCOPE( vtend         , ::vtend );

  real2d qneg("qneg",nzm,ncrms);
  real2d qpoz("poz" ,nzm,ncrms);
  int2d  nneg("nneg",nzm,ncrms);

  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qpoz(k,icrm) = 0.0;
    qneg(k,icrm) = 0.0;
    nneg(k,icrm) = 0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    t(k, j+offy_s, i+offx_s, icrm) = t(k, j+offy_s, i+offx_s, icrm) + ttend(k,icrm) * dtn;
    micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm) = 
          micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm) + qtend(k,icrm) * dtn;

    if (micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm) < 0.0) {
      yakl::atomicAdd(nneg(k,icrm),1);
      yakl::atomicAdd(qneg(k,icrm),micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm));
    } else {
      yakl::atomicAdd(qpoz(k,icrm),micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm));
    }
    dudt(na-1,k,j,i,icrm) = dudt(na-1,k,j,i,icrm) + utend(k,icrm);
    dvdt(na-1,k,j,i,icrm) = dvdt(na-1,k,j,i,icrm) + vtend(k,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real factor;
    if(nneg(k,icrm) > 0 && qpoz(k,icrm)+qneg(k,icrm) > 0.0) {
      factor =  1.0 + qneg(k,icrm)/qpoz(k,icrm);
      micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm) = 
            max(0.0,micro_field(index_water_vapor, k, j+offy_s, i+offx_s, icrm)*factor);
    }
  });
}
