
#pragma once

#include "samxx_const.h"
#include "vars.h"
#include "shear_prod2D.h"
#include "shear_prod3D.h"
#include "sat.h"

void tke_full(real5d &tke, int ind_tke, real5d &tk, int ind_tk, real5d &tkh, int ind_tkh);

