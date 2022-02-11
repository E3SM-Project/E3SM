
#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_f90.hpp"

#include "microphysics.h"
#include "micro_p3.h"

void shoc_initialize();
void shoc_proc();

enum {
  shoc_idx_qv = 0,  // total water (qv + qc)
  shoc_idx_qi,      // cloud ice amount
  shoc_idx_qr,      // rain amount
  shoc_idx_nc,      // cloud liq number
  shoc_idx_ni,      // cloud ice number
  shoc_idx_nr,      // rain number
  shoc_idx_qm,      // ice rime amount
  shoc_idx_bm,      // ice rime volume
  shoc_idx_qc,      // cloud liq amount (this is also part of total water)
  shoc_idx_tke,     // turbulent kinetic energy (tke)
  shoc_idx_esmt_u,  // explicit scalar momentum transport scheme u wind
  shoc_idx_esmt_v   // explicit scalar momentum transport scheme v wind
};
