
#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

enum {
   idx_qt = 0,  // total water (qv + qc)
   // idx_qc,      // cloud liq amount
   idx_qi,      // cloud ice amount
   idx_qr,      // rain amount
   idx_nc,      // cloud liq number
   idx_ni,      // cloud ice number
   idx_nr,      // rain number
   idx_qm,      // ice rime amount
   idx_bm       // ice rime volume
};

void micro_p3_init();
void micro_p3_proc();


