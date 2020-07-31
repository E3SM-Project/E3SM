
#pragma once

#include "samxx_const.h"
#include "vars.h"
#include "tke_full.h"
#include "diffuse_mom.h"
#include "microphysics.h"
#include "diffuse_scalar.h"

void kurant_sgs( real &cfl );

void sgs_proc();

void sgs_mom();

void sgs_scalars();

void sgs_init();

