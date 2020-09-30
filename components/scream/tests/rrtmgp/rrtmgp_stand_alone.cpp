#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/share/physics_only_grids_manager.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat.hpp"

/*
 * This will eventually contain a standalone test for the RRTMGP driver
 * As of now, it is just a shell that at least requires RRTMGP to be built
 * with the SCREAM build and test system.
 */

namespace scream {

    // Add the RRTMGP stand-alone driver test
    TEST_CASE("rrtmgp_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        /* 
         * Setup driver stuff
         */

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

        // Create the driver
        AtmosphereDriver ad;

        // Dummy timestamp
        util::TimeStamp time (0,0,0,0);

        // Initialize the driver, run the driver, cleanup
        ad.initialize(atm_comm, ad_params, time);
        ad.run(300.0);
        ad.finalize();

        // Run RRTMGP standalone codes and compare with AD run
        // Do something interesting here...
        // NOTE: these will get replaced with AD stuff that handles these
        //rrtmgp::rrtmgp_init();
        //rrtmgp::rrtmgp_main();
        //rrtmgp::rrtmgp_finalize();

        // If we got here, we were able to run the above code
        REQUIRE(true);
    }
}
