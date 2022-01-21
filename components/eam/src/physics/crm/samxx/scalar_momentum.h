#pragma once

#include "samxx_const.h"
#include "YAKL_fft.h"
#include "vars.h"

extern "C" void rfft1i(int n, real* wsave, int lensav, int ier);
extern "C" void rfft1f(int n, int inc, real* r, int lenr, real* wsave, int lensav, real* work, int lenwrk, int ier );
extern "C" void rfft1b(int n, int inc, real* r, int lenr, real* wsave, int lensav, real* work, int lenwrk, int ier );

void scalar_momentum_tend();

void scalar_momentum_pgf( real4d& scalar_wind, real4d& tend );
