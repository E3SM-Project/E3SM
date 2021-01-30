#pragma once
#include "const.h"
#include "YAKL.h"
#include "mo_gas_concentrations.h"

void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, real2d &col_dry, int ncol);


void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol);


void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol);


