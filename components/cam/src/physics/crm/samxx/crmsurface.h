
#pragma once

#include "samxx_const.h"
#include "vars.h"

void crmsurface(real1d &bflx);

YAKL_INLINE real diag_ustar(real z, real bflx, real wnd, real z0) {
  real constexpr vonk = 0.4;
  real constexpr eps = 1.0e-10;
  real constexpr am = 4.8;
  real constexpr bm = 19.3;
  real constexpr pi = 3.14159;

  real lnz = log(z/z0);
  real klnz = vonk/lnz;
  real c1 = pi/2.0 - 3.0*log(2.0);
  real ustar = wnd*klnz;

  if (bflx != 0.0) {
    for (int iterate = 0; iterate < 8; iterate++) {
      real rlmo = -bflx * vonk/(ustar*ustar*ustar + eps);
      real zeta = min(1.0,z*rlmo);
      if (zeta>0.0) {
        ustar = vonk*wnd / (lnz+am*zeta);
      } else {
        real x = sqrt(sqrt(1.0-bm*zeta));
        real psi1 = 2.0*log(1.0+x) + log(1.0+x*x) - 2.0*atan(x) + c1;
        ustar = wnd*vonk/(lnz-psi1);
      }
    }
  }
  return ustar;
}

