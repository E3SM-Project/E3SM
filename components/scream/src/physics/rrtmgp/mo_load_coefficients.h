#ifndef MO_LOAD_COEFFICIENTS_HPP
#define MO_LOAD_COEFFICIENTS_HPP

#include "cpp/rrtmgp_const.h"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include <string>

void load_and_init(GasOpticsRRTMGP &kdist, std::string filename, GasConcs const &available_gases);

#endif
