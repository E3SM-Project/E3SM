#include <iostream>
#include <cmath>
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "physics/rrtmgp/mo_garand_atmos_io.h"
#include "YAKL/Intrinsics.h"
#include "physics/rrtmgp/tests/rrtmgp_test_utils.hpp"
#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

using namespace scream;

int main (int argc, char** argv) {

    // Get filenames from command line
    if (argc != 3) {
        std::cout <<
            argv[0] << " [options] inputfile baseline\n"
            "Options:\n"
            "  (there are no options)\n";
        return 1;
    }
    std::string inputfile(argv[argc-2]);
    std::string baseline(argv[argc-1]);

    // Initialize yakl
    yakl::init();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    real2d sw_flux_up_ref;
    real2d sw_flux_dn_ref;
    real2d sw_flux_dn_dir_ref;
    real2d lw_flux_up_ref;
    real2d lw_flux_dn_ref;
    std::cout << "read_fluxes..." << std::endl;
    rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

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
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

    // Initialize the RRTMGP interface; this will read in the k-distribution
    // data that contains information about absorption coefficients for gases
    std::cout << "rrtmgp_initialize..." << std::endl;
    rrtmgp::rrtmgp_initialize(gas_concs);

    // Setup dummy all-sky problem
    real2d sfc_alb_dir;
    real2d sfc_alb_dif;
    real1d mu0;
    real2d lwp;
    real2d iwp;
    real2d rel;
    real2d rei;
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
    real2d sw_flux_dn_dir("sw_flux_dn_dir",ncol,nlay+1);
    real2d lw_flux_up ("lw_flux_up" ,ncol,nlay+1);
    real2d lw_flux_dn ("lw_flux_dn" ,ncol,nlay+1);

    // Run RRTMGP standalone codes and compare with AD run
    // Do something interesting here...
    // NOTE: these will get replaced with AD stuff that handles these
    std::cout << "rrtmgp_main..." << std::endl;
    rrtmgp::rrtmgp_main(
        p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei,
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    // Write fluxes
    std::cout << "write_fluxes..." << std::endl;
    rrtmgpTest::write_fluxes(
        baseline, 
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    // Clean up from test; this is probably not necessary, these things
    // should be deallocated when they fall out of scope, but we should be
    // good citizens and clean up our mess.
    p_lay.deallocate();
    t_lay.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();
    col_dry.deallocate();
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();
    mu0.deallocate();
    lwp.deallocate();
    iwp.deallocate();
    rel.deallocate();
    rei.deallocate();
    sw_flux_up_ref.deallocate();
    sw_flux_dn_ref.deallocate();
    sw_flux_dn_dir_ref.deallocate();
    lw_flux_up_ref.deallocate();
    lw_flux_dn_ref.deallocate();
    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dn_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();

    gas_concs.reset();
    rrtmgp::rrtmgp_finalize();
    yakl::finalize();

    return 0;
}
