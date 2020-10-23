#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "const.h"

GasOpticsRRTMGP k_dist_sw;
GasOpticsRRTMGP k_dist_lw;

// Hardcode gases for now; TODO: fix this!!!!
// Convert to string1d
const char *active_gases[8] = {
    "H2O", "CO2", "O3 ", "N2O", 
    "CO ", "CH4", "O2 ", "N2 " 
};

void rrtmgp_initialize(char const *coefficients_file_sw, char const *coefficients_file_lw) {
    // Read gas optics coefficients from file
    // Need to initialize available_gases here! The only field of the
    // available_gases type that is used int he kdist initialize is
    // available_gases%gas_name, which gives the name of each gas that would be
    // present in the ty_gas_concs object. So, we can just set this here, rather
    // than trying to fully populate the ty_gas_concs object here, which would be
    // impossible from this initialization routine because I do not thing the
    // rad_cnst objects are setup yet.
    // the other tasks!
    // TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
    //gas_names = string1d("gas_names", 8);
    GasConcs available_gases;
    load_and_init(k_dist_sw, coefficients_file_sw, available_gases);
    load_and_init(k_dist_lw, coefficients_file_lw, available_gases);
}
