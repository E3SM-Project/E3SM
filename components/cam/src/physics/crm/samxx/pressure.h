
#pragma once

#include "YAKL_fft.h"
#include "press_grad.h"
#include "press_rhs.h"
#include "samxx_const.h"
#include "vars.h"

extern "C" void fftfax_crm(int n, int *ifax, real *trigs);
extern "C" void fft991_crm(real *a, real *work, real *trigs, int *ifax, int inc, int jump, int n, int lot, int isign);

void pressure();
