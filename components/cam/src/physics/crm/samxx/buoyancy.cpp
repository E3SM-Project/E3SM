#include "buoyancy.h"

void buoyancy() {
  auto &adz   = :: adz;
  auto &dwdt  = :: dwdt;
  auto &na    = :: na;
  auto &bet   = :: bet;
  auto &tabs0 = :: tabs0;
  auto &epsv  = :: epsv;
  auto &qv    = :: qv;
  auto &qv0   = :: qv0;
  auto &qcl   = :: qcl;
  auto &qci   = :: qci;
  auto &qn0   = :: qn0;
  auto &qpl   = :: qpl;
  auto &qpi   = :: qpi;
  auto &qp0   = :: qp0;
  auto &tabs  = :: tabs;
  auto &ncrms = :: ncrms;

  if (!docolumn) {
    // for (int k=0; k<nzm-1; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kp = k+1;
      real betu, betd;
      betu = adz(k,icrm)/(adz(kp,icrm)+adz(k,icrm));
      betd = adz(kp,icrm)/(adz(kp,icrm)+adz(k,icrm));

      dwdt(na-1,kp,j,i,icrm) = 
            dwdt(na-1,kp,j,i,icrm) + 
               bet(kp,icrm)*betu*
               ( tabs0(kp,icrm)*(epsv*(qv(kp,j,i,icrm)-qv0(kp,icrm))-(qcl(kp,j,i,icrm)+qci(kp,j,i,icrm)-
                                 qn0(kp,icrm)+qpl(kp,j,i,icrm)+qpi(kp,j,i,icrm)-qp0(kp,icrm))) 
               +(tabs(kp,j,i,icrm)-tabs0(kp,icrm))*(1.0+epsv*qv0(kp,icrm)-qn0(kp,icrm)-qp0(kp,icrm)) )
               +bet(k,icrm)*betd*
               ( tabs0(k,icrm)*(epsv*(qv(k,j,i,icrm)-qv0(k,icrm))-(qcl(k,j,i,icrm)+qci(k,j,i,icrm)-
                                qn0(k,icrm)+qpl(k,j,i,icrm)+qpi(k,j,i,icrm)-qp0(k,icrm)))
               +(tabs(k,j,i,icrm)-tabs0(k,icrm))*(1.0+epsv*qv0(k,icrm)-qn0(k,icrm)-qp0(k,icrm)) );
    }); 
  }
}
