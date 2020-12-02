#include "cloud.h"

void cloud(real5d &q, int ind_q, real5d &qp, int ind_qp) {
  YAKL_SCOPE( tabs  , ::tabs );
  YAKL_SCOPE( gamaz , ::gamaz );
  YAKL_SCOPE( pres  , ::pres );
  YAKL_SCOPE( qn    , ::qn );
  YAKL_SCOPE( t     , ::t );
  YAKL_SCOPE( ncrms , ::ncrms );

  real constexpr an   = 1.0/(tbgmax-tbgmin);
  real constexpr bn   = tbgmin * an;
  real constexpr ap   = 1.0/(tprmax-tprmin);
  real constexpr bp   = tprmin * ap;
  real constexpr fac1 = fac_cond+(1.0+bp)*fac_fus;
  real constexpr fac2 = fac_fus*ap;
  real constexpr ag   = 1.0/(tgrmax-tgrmin);

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    q(ind_q,k,j+offy_s,i+offx_s,icrm)=max(0.0,q(ind_q,k,j+offy_s,i+offx_s,icrm));
    // Initial guess for temperature assuming no cloud water/ice:
    tabs(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm)-gamaz(k,icrm);
    real tabs1=(tabs(k,j,i,icrm)+fac1*qp(ind_qp,k,j+offy_s,i+offx_s,icrm))/
               (1.0+fac2*qp(ind_qp,k,j+offy_s,i+offx_s,icrm));

    real qsatt;
    real om;
    real qsatt1, qsatt2;
    // Warm cloud:
    if(tabs1 > tbgmax) {
      tabs1=tabs(k,j,i,icrm)+fac_cond*qp(ind_qp,k,j+offy_s,i+offx_s,icrm);
      qsatw_crm(tabs1,pres(k,icrm),qsatt);
    }
    // Ice cloud:
    else if(tabs1 <= tbgmin) {
      tabs1=tabs(k,j,i,icrm)+fac_sub*qp(ind_qp,k,j+offy_s,i+offx_s,icrm);
      qsati_crm(tabs1,pres(k,icrm),qsatt);
    }
    // Mixed-phase cloud:
    else {
      om = an*tabs1-bn;
       
      qsatw_crm(tabs1,pres(k,icrm),qsatt1);
      qsati_crm(tabs1,pres(k,icrm),qsatt2);
      qsatt = om*qsatt1+(1.-om)*qsatt2;
    }

    int niter;
    real dtabs, lstarn, dlstarn, omp, lstarp, dlstarp, fff, dfff, dqsat;
    //  Test if condensation is possible:
    if(q(ind_q,k,j+offy_s,i+offx_s,icrm) > qsatt) {
      niter=0;
      dtabs = 100.0;
      do {
        if(tabs1 >= tbgmax) {
          om=1.0;
          lstarn=fac_cond;
          dlstarn=0.0;
          qsatw_crm(tabs1,pres(k,icrm),qsatt);
          dtqsatw_crm(tabs1,pres(k,icrm),dqsat);
        }
        else if(tabs1 <= tbgmin) {
          om=0.0;
          lstarn=fac_sub;
          dlstarn=0.0;
          qsati_crm(tabs1,pres(k,icrm),qsatt);
          dtqsati_crm(tabs1,pres(k,icrm),dqsat);
        }
        else {
          om=an*tabs1-bn;
          lstarn=fac_cond+(1.0-om)*fac_fus;
          dlstarn=an*fac_fus;
          qsatw_crm(tabs1,pres(k,icrm),qsatt1);
          qsati_crm(tabs1,pres(k,icrm),qsatt2);
          qsatt=om*qsatt1+(1.-om)*qsatt2;
          dtqsatw_crm(tabs1,pres(k,icrm),qsatt1);
          dtqsati_crm(tabs1,pres(k,icrm),qsatt2);
          dqsat=om*qsatt1+(1.-om)*qsatt2;
        }

        if(tabs1 >= tprmax) {
          omp=1.0;
          lstarp=fac_cond;
          dlstarp=0.0;
        }
        else if(tabs1 <= tprmin) {
          omp=0.0;
          lstarp=fac_sub;
          dlstarp=0.0;
        }
        else {
          omp=ap*tabs1-bp;
          lstarp=fac_cond+(1.0-omp)*fac_fus;
          dlstarp=ap*fac_fus;
        }
        fff = tabs(k,j,i,icrm)-tabs1+lstarn*(q(ind_q,k,j+offy_s,i+offx_s,icrm)-qsatt)+
              lstarp*qp(ind_qp,k,j+offy_s,i+offx_s,icrm);
        dfff=dlstarn*(q(ind_q,k,j+offy_s,i+offx_s,icrm)-qsatt)+dlstarp*qp(ind_qp,k,j+offy_s,i+offx_s,icrm)-
             lstarn*dqsat-1.0;
        dtabs=-fff/dfff;
        niter=niter+1;
        tabs1=tabs1+dtabs;
      } while(abs(dtabs) > 0.01 && niter < 10);
      qsatt = qsatt + dqsat * dtabs;
      qn(k,j,i,icrm) = max(0.0,q(ind_q,k,j+offy_s,i+offx_s,icrm)-qsatt);
    }
    else {
      qn(k,j,i,icrm) = 0.0;
    }
    tabs(k,j,i,icrm) = tabs1;
    qp(ind_qp,k,j+offy_s,i+offx_s,icrm) = max(0.0,qp(ind_qp,k,j+offy_s,i+offx_s,icrm)); // just in case
  });

}
