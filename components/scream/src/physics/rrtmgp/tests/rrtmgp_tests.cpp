#include <iostream>
#include <cmath>
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "rrtmgp_test_utils.hpp"
#include "mo_gas_concentrations.h"
#include "mo_fluxes.h"
#include "mo_cloud_optics.h"
#include "mo_garand_atmos_io.h"
#include "Intrinsics.h"

template <class T> double arrmin(T &arr) {
    double minval = arr.myData[0];
    for (int i = 0; i<arr.totElems(); i++) {
        if (arr.myData[i] < minval) {
            minval = arr.myData[i];
        }
    }
    return minval;
}

int main(int argc, char** argv) {

    // Parse command line arguments
    if (argc != 3) {
        std::cout << "argc: " << argc << std::endl;
        std::cout << argv[0] << " <inputfile> <baseline file>\n";
        return 1;
    }
    std::string inputfile(argv[argc-2]);
    std::string baseline(argv[argc-1]);

    // Check to see that inputfiles exist
    if (!rrtmgpTest::file_exists(inputfile.c_str())) {
        std::cout << "Inputfile " << inputfile << " does not exist." << std::endl;
        return -1;
    }
    if (!rrtmgpTest::file_exists(baseline.c_str())) {
        std::cout << "Baseline " << baseline << " does not exist." << std::endl;
        return -1;
    }

    // Initialize yakl
    yakl::init();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    real2d sw_flux_up_ref;
    real2d sw_flux_dn_ref;
    real2d sw_flux_dir_ref;
    real2d lw_flux_up_ref;
    real2d lw_flux_dn_ref;
    rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

    // Get dimension sizes
    int ncol = sw_flux_up_ref.dimension[0];
    int nlev = sw_flux_up_ref.dimension[1];
    int nlay = nlev - 1;

    // Read in dummy Garand atmosphere; if this were an actual model simulation, 
    // these would be passed as inputs to the driver
    // NOTE: set ncol to size of col_flx dimension in the input file. This is so
    // that we can compare to the reference data provided in that file. Note that
    // this will copy the first column of the input data (the first profile) ncol
    // times. We will then fill some fraction of these columns with clouds for
    // the test problem.
    real2d p_lay;
    real2d t_lay;
    real2d p_lev;
    real2d t_lev;
    real2d col_dry;
    GasConcs gas_concs;
    real2d sfc_alb_dir;
    real2d sfc_alb_dif;
    real1d mu0;
    real2d lwp;
    real2d iwp;
    real2d rel;
    real2d rei;
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

    // Initialize absorption coefficients
    int ngas = gas_concs.get_num_gases();
    string1d gas_names_1d = gas_concs.get_gas_names();
    std::string gas_names[ngas];
    for (int igas = 0; igas < ngas; igas++) {
        gas_names[igas] = gas_names_1d(igas+1);
    }
    scream::rrtmgp::rrtmgp_initialize(ngas, gas_names);

    // Setup our dummy atmosphere based on the input data we read in
    rrtmgpTest::dummy_atmos(
            inputfile, ncol, p_lay, t_lay,
            sfc_alb_dir, sfc_alb_dif, mu0,
            lwp, iwp, rel, rei
        );

    // Setup flux outputs; In a real model run, the fluxes would be
    // input/outputs into the driver (persisting between calls), and
    // we would just have to setup the pointers to them in the
    // FluxesBroadband object
    real2d sw_flux_up ("sw_flux_up" ,ncol,nlay+1);
    real2d sw_flux_dn ("sw_flux_dn" ,ncol,nlay+1);
    real2d sw_flux_dir("sw_flux_dir",ncol,nlay+1);
    real2d lw_flux_up ("lw_flux_up" ,ncol,nlay+1);
    real2d lw_flux_dn ("lw_flux_dn" ,ncol,nlay+1);

    // Run RRTMGP code on dummy atmosphere
    scream::rrtmgp::rrtmgp_main(
            p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, 
            sfc_alb_dir, sfc_alb_dif, mu0,
            lwp, iwp, rel, rei,
            sw_flux_up, sw_flux_dn, sw_flux_dir,
            lw_flux_up, lw_flux_dn);

    // Check values against baseline
    rrtmgpTest::read_fluxes(
        baseline, 
        sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dir_ref,
        lw_flux_up_ref, lw_flux_dn_ref
    );
    int nerr = 0;
    if (!rrtmgpTest::all_equals(sw_flux_up_ref , sw_flux_up )) nerr++;
    if (!rrtmgpTest::all_equals(sw_flux_dn_ref , sw_flux_dn )) nerr++;
    if (!rrtmgpTest::all_equals(sw_flux_dir_ref, sw_flux_dir)) nerr++;
    if (!rrtmgpTest::all_equals(lw_flux_up_ref , lw_flux_up )) nerr++;
    if (!rrtmgpTest::all_equals(lw_flux_dn_ref , lw_flux_dn )) nerr++;

    // Clean up or else YAKL will throw errors
    scream::rrtmgp::rrtmgp_finalize();
    sw_flux_up_ref.deallocate();
    sw_flux_dn_ref.deallocate();
    sw_flux_dir_ref.deallocate();
    lw_flux_up_ref.deallocate();
    lw_flux_dn_ref.deallocate();
    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();
    p_lay.deallocate();
    t_lay.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();
    col_dry.deallocate();
    gas_concs.reset();
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();
    mu0.deallocate();
    lwp.deallocate();
    iwp.deallocate();
    rel.deallocate();
    rei.deallocate();
    yakl::finalize();

    return nerr != 0 ? 1 : 0;

}  // end of main driver code

