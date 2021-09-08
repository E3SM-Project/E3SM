
#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

enum {
   ixqv = 0,
   ixcldliq,      //cloud liquid amount index, qc
   ixcldice,      // ice index qi
   ixnumliq,      // cloud liquid number index nc
   ixnumice,      // cloud ice number index ni
   ixrain,        // rain index  qr
   ixnumrain,     // rain number index nr
   ixcldrim,      // rime index    qm
   ixrimvol       // rime volume index bm
};

void micro_p3_proc();



