
#pragma once

#include "samxx_const.h"
#include "vars.h"


YAKL_INLINE real esati_crm(real t) {

  real const a0 = 6.11147274;
  real const a1 = 0.503160820;
  real const a2 = 0.188439774e-1;
  real const a3 = 0.420895665e-3;
  real const a4 = 0.615021634e-5;
  real const a5 = 0.602588177e-7;
  real const a6 = 0.385852041e-9;
  real const a7 = 0.146898966e-11;
  real const a8 = 0.252751365e-14;

  real dtt = t-273.16;
  real esati;
  if(dtt > -80.0) {
    esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    esati = 0.01*exp(9.550426 - 5723.265/t + 3.53068*log(t) - 0.00728332*t);
  }

  return esati;
}

YAKL_INLINE real esatw_crm(real t) {

  real const a0 = 6.105851;
  real const a1 = 0.4440316;
  real const a2 = 0.1430341e-1;
  real const a3 = 0.2641412e-3;
  real const a4 = 0.2995057e-5;
  real const a5 = 0.2031998e-7;
  real const a6 = 0.6936113e-10;
  real const a7 = 0.2564861e-13;
  real const a8 = -0.3704404e-15;

  real dtt = t-273.16;

  real esatw;
  if(dtt > -80.0) {
    esatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    esatw = 2.0*0.01*exp(9.550426 - 5723.265/t + 3.53068*log(t) - 0.00728332*t);
  }
  return esatw;
}

YAKL_INLINE real dtesati_crm(real t) {

  real const a0 = 0.503223089;
  real const a1 = 0.377174432e-1;
  real const a2 = 0.126710138e-2;
  real const a3 = 0.249065913e-4;
  real const a4 = 0.312668753e-6;
  real const a5 = 0.255653718e-8;
  real const a6 = 0.132073448e-10;
  real const a7 = 0.390204672e-13;
  real const a8 = 0.497275778e-16;

  real dtt = t-273.16;
  real dtesati;
  if(dtt > -80.0) {
    dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    dtesati= esati_crm(t+1.0)-esati_crm(t);
  }

  return dtesati;
}

YAKL_INLINE real dtesatw_crm(real t) {

  real const a0 = 0.443956472;
  real const a1 = 0.285976452e-1;
  real const a2 = 0.794747212e-3;
  real const a3 = 0.121167162e-4;
  real const a4 = 0.103167413e-6;
  real const a5 = 0.385208005e-9;
  real const a6 = -0.604119582e-12;
  real const a7 = -0.792933209e-14;
  real const a8 = -0.599634321e-17;

  real dtt = t-273.16;
  real dtesatw;
  if(dtt > -80.0) {
    dtesatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    dtesatw = esatw_crm(t+1.0)-esatw_crm(t);
  }

  return dtesatw;
}

YAKL_INLINE void qsati_crm(real t, real p, real &qsati) {

  real esati;
  esati = esati_crm(t);
  qsati = 0.622*esati/max(esati,p-esati);
}

YAKL_INLINE void qsatw_crm(real t, real p, real &qsatw) {

  real esatw;
  esatw = esatw_crm(t);
  qsatw = 0.622*esatw/max(esatw,p-esatw);
}
YAKL_INLINE void dtqsati_crm(real t, real p, real &dtqsati) {

  dtqsati = 0.622*dtesati_crm(t)/p;

}

YAKL_INLINE void dtqsatw_crm(real t, real p, real &dtqsatw) {

  dtqsatw = 0.622*dtesatw_crm(t)/p;

}

