
#pragma once

#include "samxx_const.h"
#include "vars.h"
#include "cloud.h"
#include "precip_init.h"
#include "precip_proc.h"

YAKL_INLINE real term_vel_qp(int icrm, int i, int j, int k, real qploc, real vrain, real vsnow, 
                             real vgrau, real crain, real csnow, real cgrau, real rho,
                             real tabs, real a_pr, real a_gr) {
  real term_vel = 0.0;
  if(qploc > qp_threshold) {
    real omp = max(0.0,min(1.0,(tabs-tprmin)*a_pr));
    if(omp == 1.0) {
      term_vel = vrain*pow(rho*qploc,crain);
    }
    else if(omp == 0.0) {
      real omg = max(0.0,min(1.0,(tabs-tgrmin)*a_gr));
      real qgg=omg*qploc;
      real qss=qploc-qgg;
      term_vel = (omg*vgrau*pow(rho*qgg,cgrau) + (1.0-omg)*vsnow*pow(rho*qss,csnow));
    }
    else {
      real omg = max(0.0,min(1.0,(tabs-tgrmin)*a_gr));
      real qrr=omp*qploc;
      real qss=qploc-qrr;
      real qgg=omg*qss;
      qss=qss-qgg;
      term_vel = (omp*vrain*pow(rho*qrr,crain) + (1.0-omp)*(omg*vgrau*pow(rho*qgg,cgrau) + (1.0-omg)*vsnow*pow(rho*qss,csnow)));
    }
  }
  return term_vel;

}

void precip_fall(int hydro_type, real4d &omega);

void micro_precip_fall();

void micro_flux();

void micro_diagnose();

void micro_proc();

void micro_init();

YAKL_INLINE real pp(real y) {
  return max(0.0,y);
}

YAKL_INLINE real pn(real y) {
  return -min(0.0,y);
}


