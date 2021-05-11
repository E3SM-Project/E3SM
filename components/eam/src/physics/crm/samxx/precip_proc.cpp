
#include "precip_proc.h"

void precip_proc(real5d &q, int ind_q, real5d &qp, int ind_qp) {
  auto &tabs          = :: tabs;
  auto &a_bg          = :: a_bg;
  auto &a_pr          = :: a_pr;
  auto &a_gr          = :: a_gr;
  auto &coefice       = :: coefice;
  auto &accrrc        = :: accrrc;
  auto &accrsc        = :: accrsc;
  auto &accrsi        = :: accrsi;
  auto &accrgc        = :: accrgc;
  auto &accrgi        = :: accrgi;
  auto &dtn           = :: dtn;
  auto &pres          = :: pres;
  auto &evapr1        = :: evapr1;
  auto &evapr2        = :: evapr2;
  auto &evaps1        = :: evaps1;
  auto &evaps2        = :: evaps2;
  auto &evapg1        = :: evapg1;
  auto &evapg2        = :: evapg2;
  auto &qpsrc         = :: qpsrc;
  auto &qpevp         = :: qpevp;
  auto &qn            = :: qn;
  auto &ncrms         = :: ncrms;

  real powr1 = (3.0 + b_rain) / 4.0;
  real powr2 = (5.0 + b_rain) / 8.0;
  real pows1 = (3.0 + b_snow) / 4.0;
  real pows2 = (5.0 + b_snow) / 8.0;
  real powg1 = (3.0 + b_grau) / 4.0;
  real powg2 = (5.0 + b_grau) / 8.0;

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qpsrc(k,icrm)=0.0;
    qpevp(k,icrm)=0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    //-------     Autoconversion/accretion
    real omn, omp, omg, qcc, qii, autor, autos, accrr, qrr, accrcs, accris,
         qss, accrcg, accrig, tmp, qgg, dq, qsatt, qsat;

    if (qn(k,j,i,icrm)+qp(ind_qp,k,j+offy_s,i+offx_s,icrm) > 0.0) {
      omn = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tbgmin)*a_bg));
      omp = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tprmin)*a_pr));
      omg = max(0.0,min(1.0,(tabs(k,j,i,icrm)-tgrmin)*a_gr));

      if (qn(k,j,i,icrm) > 0.0) {
        qcc = qn(k,j,i,icrm) * omn;
        qii = qn(k,j,i,icrm) * (1.0-omn);

        if (qcc > qcw0) {
          autor = alphaelq;
        } else {
          autor = 0.0;
        }

        if (qii > qci0) {
          autos = betaelq*coefice(k,icrm);
        } else {
          autos = 0.0;
        }

        accrr = 0.0;
        if (omp > 0.001) {
          qrr = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * omp;
          accrr = accrrc(k,icrm) * pow(qrr, powr1);
        }

        accrcs = 0.0;
        accris = 0.0;

        if (omp < 0.999 && omg < 0.999) {
          qss = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * (1.0-omp)*(1.0-omg);
          tmp = pow(qss, pows1);
          accrcs = accrsc(k,icrm) * tmp;
          accris = accrsi(k,icrm) * tmp;
        }
        accrcg = 0.0;
        accrig = 0.0;
        if (omp < 0.999 && omg > 0.001) {
          qgg = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * (1.0-omp)*omg;
          tmp = pow(qgg, powg1);
          accrcg = accrgc(k,icrm) * tmp;
          accrig = accrgi(k,icrm) * tmp;
        }
        qcc = (qcc+dtn*autor*qcw0)/(1.0+dtn*(accrr+accrcs+accrcg+autor));
        qii = (qii+dtn*autos*qci0)/(1.0+dtn*(accris+accrig+autos));
        dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+(accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0));
        dq = min(dq,qn(k,j,i,icrm));
        qp(ind_qp,k,j+offy_s,i+offx_s,icrm) = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) + dq;
        q(ind_q,k,j+offy_s,i+offx_s,icrm) = q(ind_q,k,j+offy_s,i+offx_s,icrm) - dq;
        qn(k,j,i,icrm) = qn(k,j,i,icrm) - dq;
        yakl::atomicAdd(qpsrc(k,icrm),dq);

      } else if(qp(ind_qp,k,j+offy_s,i+offx_s,icrm) > qp_threshold && qn(k,j,i,icrm) == 0.0) {

        qsatt = 0.0;
        if(omn > 0.001) {
          qsatw_crm(tabs(k,j,i,icrm),pres(k,icrm),qsat);
          qsatt = qsatt + omn*qsat;
        }
        if(omn < 0.999) {
          qsati_crm(tabs(k,j,i,icrm),pres(k,icrm),qsat);
          qsatt = qsatt + (1.-omn)*qsat;
        }
        dq = 0.0;
        if(omp > 0.001) {
          qrr = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * omp;
          dq = dq + evapr1(k,icrm)*sqrt(qrr) + evapr2(k,icrm)*pow(qrr,powr2);
        }
        if(omp < 0.999 && omg < 0.999) {
          qss = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * (1.0-omp)*(1.0-omg);
          dq = dq + evaps1(k,icrm)*sqrt(qss) + evaps2(k,icrm)*pow(qss,pows2);
        }
        if(omp < 0.999 && omg > 0.001) {
          qgg = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) * (1.0-omp)*omg;
          dq = dq + evapg1(k,icrm)*sqrt(qgg) + evapg2(k,icrm)*pow(qgg,powg2);
        }
        dq = dq * dtn * (q(ind_q,k,j+offy_s,i+offx_s,icrm) /qsatt-1.0);
        dq = max(-0.5*qp(ind_qp,k,j+offy_s,i+offx_s,icrm),dq);
        qp(ind_qp,k,j+offy_s,i+offx_s,icrm) = qp(ind_qp,k,j+offy_s,i+offx_s,icrm) + dq;
        q(ind_q,k,j+offy_s,i+offx_s,icrm) = q(ind_q,k,j+offy_s,i+offx_s,icrm) - dq;
        yakl::atomicAdd(qpevp(k,icrm),dq);

      } else { 

        q(ind_q,k,j+offy_s,i+offx_s,icrm) = q(ind_q,k,j+offy_s,i+offx_s,icrm) + qp(ind_qp,k,j+offy_s,i+offx_s,icrm);
        yakl::atomicAdd(qpevp(k,icrm),-qp(ind_qp,k,j+offy_s,i+offx_s,icrm));
        qp(ind_qp,k,j+offy_s,i+offx_s,icrm) = 0.0;

      }

    }

    dq = qp(ind_qp,k,j+offy_s,i+offx_s,icrm);
    qp(ind_qp,k,j+offy_s,i+offx_s,icrm)=max(0.0,qp(ind_qp,k,j+offy_s,i+offx_s,icrm));
    q(ind_q,k,j+offy_s,i+offx_s,icrm) = q(ind_q,k,j+offy_s,i+offx_s,icrm) + (dq-qp(ind_qp,k,j+offy_s,i+offx_s,icrm));

  });
}


