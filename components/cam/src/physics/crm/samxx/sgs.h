
#pragma once

#include "diffuse_mom.h"
#include "diffuse_scalar.h"
#include "microphysics.h"
#include "samxx_const.h"
#include "tke_full.h"
#include "vars.h"

void kurant_sgs(real &cfl);

void sgs_proc();

void sgs_mom();

void sgs_scalars();

void sgs_init();
