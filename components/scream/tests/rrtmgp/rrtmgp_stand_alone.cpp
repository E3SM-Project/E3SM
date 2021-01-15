#include <iostream>
#include <cmath>
#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/share/physics_only_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "Intrinsics.h"
#include "rrtmgp_test_utils.hpp"

/*
 * Run standalone test problem for RRTMGP. Two tests are run, one that uses
 * just the RRTMGP interface, one that runs the test through the SCREAM AD to
 * test the coupling with SCREAM. These should be produce identical results,
 * and the results from both tests are each compared with the reference fluxes
 * saved in the input file.
 */

// Input file that contains example atmosphere and reference fluxes
std::string inputfile = "./data/rrtmgp-allsky.nc";

namespace scream {

    // Add the RRTMGP stand-alone driver test
    TEST_CASE("rrtmgp_stand_alone", "") {

         // Setup for standalone (dummy) problem

        // Initialize yakl
        yakl::init();

        // Get reference fluxes from input file; do this here so we can get ncol dimension
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

        // Get dimension sizes
        int ncol = sw_flux_up_ref.dimension[0];
        int nlev = sw_flux_up_ref.dimension[1];
        int nlay = nlev - 1;

        // Read atmosphere profile
        real2d p_lay;
        real2d t_lay;
        real2d p_lev;
        real2d t_lev;
        real2d col_dry;
        GasConcs gas_concs;
        read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

        // Initialize the RRTMGP interface; this will read in the k-distribution
        // data that contains information about absorption coefficients for gases
        int ngas = gas_concs.get_num_gases();
        string1d gas_names_1d = gas_concs.get_gas_names();
        std::string gas_names[ngas];
        for (int igas = 0; igas < ngas; igas++) {
            gas_names[igas] = gas_names_1d(igas+1);
        }
        rrtmgp::rrtmgp_initialize(ngas, gas_names);

        // Setup a dummy all-sky atmosphere; if this were an actual model simulation,
        // these would be passed as inputs to the driver
        // NOTE: set ncol to size of col_flx dimension in the input file. This is so
        // that we can compare to the reference data provided in that file. Note that
        // this will copy the first column of the input data (the first profile) ncol
        // times. We will then fill some fraction of these columns with clouds for
        // the test problem.
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
        rrtmgp::rrtmgp_main(
            p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry,
            sfc_alb_dir, sfc_alb_dif, mu0,
            lwp, iwp, rel, rei,
            sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
            lw_flux_up, lw_flux_dn
        );

        // Check values
        REQUIRE(rrtmgpTest::all_equals(sw_flux_up_ref    , sw_flux_up    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_ref    , sw_flux_dn    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_dir_ref, sw_flux_dn_dir));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_up_ref    , lw_flux_up    ));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_dn_ref    , lw_flux_dn    ));

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

        // Make sure we exit cleanly
        REQUIRE(true);
    }
        
    /* 
     * Run standalone test through SCREAM driver this time
     */
    TEST_CASE("rrtmgp_scream_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        // Initialize yakl
        yakl::init();

        // Read reference fluxes from input file
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

        // Get dimension sizes
        int ncol = sw_flux_up_ref.dimension[0];
        int nlev = sw_flux_up_ref.dimension[1];
        int nlay = nlev - 1;

        // Load ad parameter list
        std::string fname = "input.yaml";
        ekat::ParameterList ad_params("Atmosphere Driver");
        REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

        // Create a MPI communicator
        ekat::Comm atm_comm (MPI_COMM_WORLD);

        // Need to register products in the factory *before* we create any atm process or grids manager.,
        auto& proc_factory = AtmosphereProcessFactory::instance();
        auto& gm_factory = GridsManagerFactory::instance();
        proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);
        gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);

        // Create the grids manager
        auto& gm_params = ad_params.sublist("Grids Manager");
        const std::string& gm_type = gm_params.get<std::string>("Type");
        auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

        // Create the driver
        AtmosphereDriver ad;

        // Dummy timestamp
        util::TimeStamp time (0,0,0,0);

        // Initialize the driver, run the driver, cleanup
        ad.initialize(atm_comm, ad_params, time);
        ad.run(300.0);

        // Check values; need to get fluxes from field manager first
        // The AD should have called RRTMGP to calculate these values in the ad.run() call
        auto& field_repo = ad.get_field_repo();
        auto& d_sw_flux_up = field_repo.get_field("sw_flux_up", "Physics").get_view();
        auto& d_sw_flux_dn = field_repo.get_field("sw_flux_dn", "Physics").get_view();
        auto& d_sw_flux_dn_dir = field_repo.get_field("sw_flux_dn_dir", "Physics").get_view();
        auto& d_lw_flux_up = field_repo.get_field("lw_flux_up", "Physics").get_view();
        auto& d_lw_flux_dn = field_repo.get_field("lw_flux_dn", "Physics").get_view();
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_up_test("sw_flux_up_test", d_sw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn_test("sw_flux_dn_test", d_sw_flux_dn.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn_dir_test("sw_flux_dn_dir_test", d_sw_flux_dn_dir.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_up_test("lw_flux_up_test", d_lw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_dn_test("lw_flux_dn_test", d_lw_flux_dn.data(), ncol, nlay+1);

        // Make sure fluxes from field manager that were calculated in AD call of RRTMGP match reference fluxes from input file
        REQUIRE(rrtmgpTest::all_equals(sw_flux_up_ref    , sw_flux_up_test  ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_ref    , sw_flux_dn_test    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_dir_ref, sw_flux_dn_dir_test));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_up_ref    , lw_flux_up_test    ));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_dn_ref    , lw_flux_dn_test    ));

        // Clean up after ourselves
        sw_flux_up_ref.deallocate();
        sw_flux_dn_ref.deallocate();
        sw_flux_dn_dir_ref.deallocate();
        lw_flux_up_ref.deallocate();
        lw_flux_dn_ref.deallocate();
        sw_flux_up_test.deallocate();
        sw_flux_dn_test.deallocate();
        sw_flux_dn_dir_test.deallocate();
        lw_flux_up_test.deallocate();
        lw_flux_dn_test.deallocate();
        ad.finalize();
        yakl::finalize();

        // If we got this far, we were able to run the code through the AD
        REQUIRE(true);
    }
}
