
#pragma once

#include "samxx_const.h"
#include "samxx_utils.h"
#include "vars.h"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_f90.hpp"

void shoc_initialize();
void shoc_proc();
template <typename ArrayT>
void shoc_update_precipitation(const ArrayT& qtracers);

