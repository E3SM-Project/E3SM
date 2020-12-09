#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "const.h"

// Prototypes
extern "C" void add_gas_name(char const *gas_name);
extern "C" void convert_gas_names(string1d &gas_names);
extern "C" void vect_to_string1d(std::vector<std::string> vect, string1d strarr);
extern "C" int get_nbands_sw();
extern "C" int get_nbands_lw();
extern "C" void rrtmgpxx_finalize();

GasOpticsRRTMGP k_dist_sw;
GasOpticsRRTMGP k_dist_lw;

// Vector of strings to hold active gas names. 
// These need to be added at runtime, one by one,
// via the add_gas_name function.
std::vector<std::string> gas_names_vect;

extern "C" void rrtmgpxx_initialize_cpp(char const *coefficients_file_sw, char const *coefficients_file_lw) {
    // First, make sure yakl has been initialized
    if (!yakl::isInitialized()) {
        yakl::init();
    }

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
    //
    // Let us cheat for a moment and hard-code the gases.
    // TODO: fix this!
    string1d gas_names("gas_names", gas_names_vect.size());
    convert_gas_names(gas_names);
    GasConcs available_gases;
    available_gases.init(gas_names, 1, 1);
    load_and_init(k_dist_sw, coefficients_file_sw, available_gases);
    load_and_init(k_dist_lw, coefficients_file_lw, available_gases);
}

extern "C" void rrtmgpxx_finalize() {
    k_dist_sw.finalize();
    k_dist_lw.finalize();
    yakl::finalize();
}

extern "C" void add_gas_name(char const *gas_name) {
    gas_names_vect.push_back(std::string(gas_name));
}

extern "C" void convert_gas_names(string1d &gas_names) {
    int ngas = gas_names_vect.size();
    if (ngas == 0) {
        throw "No active gases; are you sure you initialized gas_names_vect?";
    }
    for (int i = 1; i <= ngas; i++) {
        gas_names(i) = gas_names_vect[i-1];
    }
}

extern "C" void vect_to_string1d(std::vector<std::string> vect, string1d strarr) {
    int n = vect.size();
    for (int i = 0; i < n; i++) {
        strarr(i+1) = vect[i];
    }
}

extern "C" int get_nbands_sw() {
    return k_dist_sw.get_nband();
}

extern "C" int get_nbands_lw() {
    return k_dist_lw.get_nband();
}
