#ifndef MO_GARAND_ATMOS_IO_HPP
#define MO_GARAND_ATMOS_IO_HPP

#include "cpp/rrtmgp_const.h"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "YAKL.h"

void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, int ncol);


void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol);


void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol);

#endif
