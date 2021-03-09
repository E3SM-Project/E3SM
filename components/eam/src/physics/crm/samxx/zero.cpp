
#include "zero.h"

void zero() {
  auto &dudt   = :: dudt;
  auto &dvdt   = :: dvdt;
  auto &dwdt   = :: dwdt;
  auto &misc   = :: misc;
  auto &na     = :: na;
  auto &ncrms  = :: ncrms;

  // for (int k=0; k<nz; k++) {
  //   for (int j=0; j<nyp1; j++) {
  //     for (int i=0; i<nxp1; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nz,nyp1,nxp1,ncrms) , YAKL_LAMBDA (int k, int j , int i, int icrm) {
    if(i<nxp1 && j<ny && k<nzm){ dudt(na-1,k,j,i,icrm) = 0.0; }
    if(i<nx && j<nyp1 && k<nzm){ dvdt(na-1,k,j,i,icrm) = 0.0; }
    if(i<nx && j<ny && k<nz){ dwdt(na-1,k,j,i,icrm) = 0.0; }
    if(i<nx && j<ny && k<nz){ misc(k,j,i,icrm) = 0.0; }
  }); 
}


