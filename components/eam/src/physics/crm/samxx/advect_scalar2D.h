
#pragma once

#include "samxx_const.h"
#include "vars.h"

void advect_scalar2D(real4d &f, real2d &flux);

void advect_scalar2D(real5d &f, int ind_f, real2d &flux);

void advect_scalar2D(real5d &f, int ind_f, real3d &flux, int ind_flux);

YAKL_INLINE real andiff2(real x1, real x2, real a, real b) {
  return (abs(a)-a*a*b)*0.5*(x2-x1);
}

YAKL_INLINE real across2(real x1, real a1, real a2) {
  return 0.03125*a1*a2*x1;
}

YAKL_INLINE real pp2(real y) {
  return max(0.0,y);
}

YAKL_INLINE real pn2(real y) {
  return -min(0.0,y);
}

