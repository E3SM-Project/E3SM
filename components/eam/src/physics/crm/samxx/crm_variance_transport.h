#pragma once

#include "samxx_const.h"
#include "vars.h"
#include "YAKL_fft.h"

void VT_filter(real4d &f_in, real4d &f_out);
void VT_diagnose();
void VT_forcing();