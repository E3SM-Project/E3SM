#pragma once

#include "samxx_const.h"
#include "YAKL_fft.h"
#include "vars.h"

extern "C" void rfft1i(int n, real* wsave, int lensav, int ier);
extern "C" void rfft1f(int n, int inc, real* r, int lenr, real* wsave, int lensav, real* work, int lenwrk, int ier );
extern "C" void rfft1b(int n, int inc, real* r, int lenr, real* wsave, int lensav, real* work, int lenwrk, int ier );

void esmt_fft_forward(real2d& arr_in, real1d& k_out, real2d& arr_out);

void esmt_fft_backward(real2d& arr_in, real2d& arr_out);

void scalar_momentum_tend();

void scalar_momentum_pgf( int icrm, real3d& u_s, real3d& tend );
