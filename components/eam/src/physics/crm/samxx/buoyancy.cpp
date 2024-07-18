#include "buoyancy.h"

void buoyancy() {
  YAKL_SCOPE( adz               , :: adz);
  YAKL_SCOPE( dwdt              , :: dwdt);
  YAKL_SCOPE( na                , :: na);
  YAKL_SCOPE( bet               , :: bet);
  YAKL_SCOPE( tabs0             , :: tabs0);
  YAKL_SCOPE( epsv              , :: epsv);
  YAKL_SCOPE( qv                , :: qv);
  YAKL_SCOPE( qv0               , :: qv0);
  YAKL_SCOPE( qcl               , :: qcl);
  YAKL_SCOPE( qci               , :: qci);
  YAKL_SCOPE( bou0              , :: bou0);
  YAKL_SCOPE( qn0               , :: qn0);
  YAKL_SCOPE( qc0               , :: qc0);
  YAKL_SCOPE( qi0               , :: qi0);
  YAKL_SCOPE( qpl               , :: qpl);
  YAKL_SCOPE( qpi               , :: qpi);
  YAKL_SCOPE( qp0               , :: qp0);
  YAKL_SCOPE( tabs              , :: tabs);
  YAKL_SCOPE( rho               , :: rho);
  YAKL_SCOPE( ncrms             , :: ncrms);
  // YAKL_SCOPE( crm_output_bou_ls , :: crm_output_bou_ls );
  // YAKL_SCOPE( crm_output_bou    , :: crm_output_bou );
  YAKL_SCOPE( crm_output_tkeqc  , :: crm_output_tkeqc );
  YAKL_SCOPE( crm_output_tkeqt  , :: crm_output_tkeqt );
  YAKL_SCOPE( crm_output_tkeb   , :: crm_output_tkeb );

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
      // crm_output_bou(k,j,i,icrm) = tabs(k,j,i,icrm)*(1+0.61*qv(k,j,i,icrm)-(qcl(k,j,i,icrm)+qci(k,j,i,icrm)+qpl(k,j,i,icrm)+qpi(k,j,i,icrm)));
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

    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) { 
      real tmp_bou = tabs(k,j,i,icrm)*(1+0.61*qv(k,j,i,icrm)-(qcl(k,j,i,icrm)+qci(k,j,i,icrm)+qpl(k,j,i,icrm)+qpi(k,j,i,icrm)));
      yakl::atomicAdd(bou0(k,icrm) , tmp_bou);
      yakl::atomicAdd(qc0(k,icrm) , qcl(k,j,i,icrm));
      yakl::atomicAdd(qi0(k,icrm) , qci(k,j,i,icrm));
    }); 

    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
        int l = plev-(k+1);
        int kp = k+1;
        real u2z = 0.0;
        real v2z = 0.0;
        real w2z = 0.0;
        real wqcz = 0.0;
        real wqtz = 0.0;
        real wb   = 0.0;
        for (int j=0; j<ny; j++) {
          for (int i=0; i<nx; i++) {
            real tmp_qt = qcl(k,j,i,icrm) + qci(k,j,i,icrm) + qpl(k,j,i,icrm) + qpi(k,j,i,icrm);
            real tmp_bou = tabs(k,j,i,icrm) * ( 1 + 0.61*qv(k,j,i,icrm) - tmp_qt );
            real tmp_qt_kp = qcl(kp,j,i,icrm) + qci(kp,j,i,icrm) + qpl(kp,j,i,icrm) + qpi(kp,j,i,icrm);
            real tmp_bou_kp = tabs(kp,j,i,icrm) * ( 1 + 0.61*qv(kp,j,i,icrm) - tmp_qt_kp );
            real tmp3 = dwdt(na-1,kp,j,i,icrm);
            real tmp4 = dwdt(na-1,k,j,i,icrm);
            real tmp5 = qcl(kp,j,i,icrm)-qc0(kp,icrm);
            real tmp6 = qcl(k,j,i,icrm)-qc0(k,icrm);
            real tmp7 = qci(kp,j,i,icrm)-qi0(kp,icrm);
            real tmp8 = qci(k,j,i,icrm)-qi0(k,icrm);
            real tmp9 = tmp5+tmp7;
            real tmp10 = tmp6+tmp8;
            // real tmp11 = tmp_bou_kp-bou0(kp,icrm);
            // real tmp12 = tmp_bou   -bou0(k,icrm);
            real tmp11 = 9.81 * ( tmp_bou_kp - bou0(kp,icrm) ) / tmp_bou_kp;
            real tmp12 = 9.81 * ( tmp_bou    - bou0(k ,icrm) ) / tmp_bou   ; 
            w2z = w2z+0.5*(tmp3*tmp3+tmp4*tmp4);
            wqcz = wqcz+0.5*(tmp3*tmp5+tmp4*tmp6);
            wqtz = wqtz+0.5*(tmp3*tmp9+tmp4*tmp10);
            wb   = wb+0.5*(tmp3*tmp11+tmp4*tmp12);
          }
        }
        crm_output_tkeqc(l,icrm)= crm_output_tkeqc(l,icrm) + rho(k,icrm)*0.5*wqcz;
        crm_output_tkeqt(l,icrm)= crm_output_tkeqt(l,icrm) + rho(k,icrm)*0.5*wqtz;
        crm_output_tkeb (l,icrm)= crm_output_tkeb(l,icrm)  + rho(k,icrm)*0.5*wb;
    }); 



  }

}
