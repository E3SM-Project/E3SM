#include <catch2/catch.hpp>
#include "share/atm_process/atmosphere_process.hpp"
#include "share/scream_pack.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "control/atmosphere_driver.hpp"

/*
 * This will eventually contain a standalone test for the RRTMGP driver 
 * As of now, it is just a shell that at least requires RRTMGP to be built
 * with the SCREAM build and test system. 
 */

#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
// #include "physics/rrtmgp/rrtmgp_functions_f90.hpp"

namespace scream {
    // A dummy physics grids for this test //
    class DummyPhysicsGrid : public SEGrid {
        public: DummyPhysicsGrid (const int num_cols) : SEGrid("Physics",GridType::SE_NodeBased,num_cols) {
            // Nothing to do here
        }
        ~DummyPhysicsGrid () = default;
    };

    // Add the RRTMGP stand-alone driver test
    TEST_CASE("rrtmgp_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        constexpr int num_iter = 20;
        constexpr int num_cols = 10;

        /* 
         * Setup driver stuff
         */

        // Setup parameter list for inputs to the radiation interface
        ParameterList ad_params("Atmosphere Driver");

        // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
        // which rely on factory for process creation. The initialize method of the AD does that.
        // While we're at it, check that the case insensitive key of the factory works.
        auto& proc_factory = AtmosphereProcessFactory::instance();
        proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);

        // Create a MPI communicator
        Comm atm_comm (MPI_COMM_WORLD);

        // Create the driver
        AtmosphereDriver ad;

        // Dummy timestamp
        util::TimeStamp time (0,0,0);

        // Initialize the driver
        ad.initialize(atm_comm, ad_params, time);

        // Do something interesting here...
        // NOTE: these will get replaced with AD stuff that handles these
        rrtmgp::rrtmgp_init();
        rrtmgp::rrtmgp_main();
        rrtmgp::rrtmgp_finalize();

        // If we got here, we were able to run the above code
        REQUIRE(true);
    }
}
