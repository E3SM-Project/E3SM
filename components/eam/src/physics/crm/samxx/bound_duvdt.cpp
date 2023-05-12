#include "bound_duvdt.h"

void bound_duvdt() {
  YAKL_SCOPE( dudt   , ::dudt);
  YAKL_SCOPE( dvdt   , ::dvdt);
  YAKL_SCOPE( na     , ::na);  
  YAKL_SCOPE( ncrms  , ::ncrms);

  // for (int k=0; k<nzm; k++) {
  //  for (int j=0; j<ny; j++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    dudt(na-1,k,j,nx,icrm) = dudt(na-1,k,j,0,icrm);
  });

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dvdt(na-1,k,ny,i,icrm) = dvdt(na-1,k,0,i,icrm);
    });
  }
}
