
#pragma once

#include "samxx_const.h"
#include "vars.h"
#include "advect_scalar2D.h"
#include "advect_scalar3D.h"

void advect_scalar(real4d &f, real2d &fadv, real2d &flux);

void advect_scalar(real5d &f, int ind_f, real2d &fadv, real2d &flux);

void advect_scalar(real5d &f, int ind_f, real3d &fadv, int ind_fadv, real3d &flux, int ind_flux);

