
#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"
// #include "p3_lookup.hpp"

void micro_p3_init();
void micro_p3_proc();
void micro_p3_diagnose();


using namespace scream;
using namespace scream::p3;
using P3F = Functions<Real, DefaultDevice>;

// Table dimension constants from scream/src/physics/share/physics_constants.hpp
const int VTABLE_DIM0    = 300;
const int VTABLE_DIM1    = 10;
const int MU_R_TABLE_DIM = 150;

// ice lookup table constants from scream/src/physics/p3/p3_functions.hpp
const int densize             = 5;
const int rimsize             = 4;
const int isize               = 50;
const int ice_table_size      = 12; // number of quantities used from lookup table
const int rcollsize           = 30;
const int collect_table_size  = 2;  // number of ice-rain collection  quantities used from lookup table
const int iparam              = 3;
const int dnusize             = 16;

extern real4d ice_table_save;      // ice lookup table values
extern real5d collect_table_save;  // ice lookup table values for ice-rain collision/collection
extern real2d vn_table_save;       // lookup table values for rain number
extern real2d vm_table_save;       // lookup table values for mass-weighted fallspeeds
extern real2d revap_table_save;    // lookup table values for ventilation parameters
extern real1d mu_r_table_save;     // lookup table values for rain shape parameter
extern real1d dnu_table_save;      // droplet spectral shape parameter for mass spectra

extern "C" void initialize_p3_lookup();
extern "C" void finalize_p3_lookup();
