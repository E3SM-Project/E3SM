
#pragma once

#include "samxx_const.h"
#include "vars.h"

YAKL_INLINE void z0_est(real z, real bflx, real wnd, real ustar, real &z0) {
  real vonk = 0.4;
  real eps = 1.0e-10;
  real am = 4.8;
  real bm = 19.3;
  real c1 = 3.14159/2.0 - 3.0*log(2.0);
  real rlmo = -bflx*vonk/(ustar*ustar*ustar+eps);
  real zeta = min(1.0,z*rlmo);

  real x;
  real psi1;
  if(zeta >= 0.0) {
    psi1 = -am*zeta;
  }
  else {
    x = sqrt(sqrt(1.0-bm*zeta));
    psi1 = 2.0*log(1.0+x) + log(1.0+x*x) -2.0*atan(x) + c1;
  }

  real lnz = max(0.0, vonk*wnd/(ustar+eps) +psi1);

  z0 = z*exp(-lnz);
}

