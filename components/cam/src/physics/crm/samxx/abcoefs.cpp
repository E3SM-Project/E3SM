
#include "abcoefs.h"

// Compute the coefficients for the Adams-Bashforth scheme
void abcoefs() {
  if (nstep >= 3) {
    realHost1d dt3Host("dt3Host",3);
    dt3.deep_copy_to(dt3Host);
    yakl::fence();
    real alpha = dt3Host(nb-1) / dt3Host(na-1);
    real beta  = dt3Host(nc-1) / dt3Host(na-1);
    ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta);
    bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha);
    at = 1. - bt - ct;
  } else if (nstep >= 2) {
    at = 3./2.;
    bt = -1./2.;
    ct = 0.;
  } else {
    at = 1.;
    bt = 0.;
    ct = 0.;
  }
}

