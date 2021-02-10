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
#include "ekat/util/ekat_test_utils.hpp"
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "Intrinsics.h"
#include "rrtmgp_test_utils.hpp"

/*
 * Run standalone test problem for RRTMGP and compare with baseline
 */

// Input file that contains example atmosphere and reference fluxes

namespace scream {
       
    /* 
     * Run standalone test through SCREAM driver this time
     */
    TEST_CASE("rrtmgp_scream_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        // Get baseline name (needs to be passed as an arg)
        //std::string inputfile = "./data/rrtmgp-allsky.nc";
        std::string inputfile = ekat::TestSession::get().params.at("rrtmgp_inputfile");
        std::string baseline = ekat::TestSession::get().params.at("rrtmgp_baseline");

        // Check if files exists
        REQUIRE(rrtmgpTest::file_exists(inputfile.c_str()));
        REQUIRE(rrtmgpTest::file_exists(baseline.c_str()));

        // Initialize yakl
        yakl::init();

        // Read reference fluxes from input file
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        // TODO: this should read BASELINE instead, but those need to exist
        // before cmake is run?
        rrtmgpTest::read_fluxes(baseline, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

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
