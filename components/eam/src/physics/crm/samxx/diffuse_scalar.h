
#pragma once

#include "samxx_const.h"
#include "vars.h"

#include "diffuse_scalar2D.h"
#include "diffuse_scalar3D.h"

void diffuse_scalar(real5d &tkh, int ind_tkh, real4d &f, real3d &fluxb, real3d &fluxt,
                    real2d &fdiff, real2d &flux);
void diffuse_scalar(real5d &tkh, int ind_tkh, real5d &f, int ind_f, real3d &fluxb,
                    real3d &fluxt, real2d &fdiff, real2d &flux);
void diffuse_scalar(real5d &tkh, int ind_tkh, real5d &f, int ind_f, real4d &fluxb, int ind_fluxb,
                    real4d &fluxt, int ind_fluxt, real3d &fdiff, int ind_fdiff, real3d &flux, int ind_flux);

