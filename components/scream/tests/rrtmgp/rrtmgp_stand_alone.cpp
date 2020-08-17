#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"

#include "ekat/ekat.hpp"

/*
 * This will eventually contain a standalone test for the RRTMGP driver
 * As of now, it is just a shell that at least requires RRTMGP to be built
 * with the SCREAM build and test system.
 */

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
        ekat::ParameterList ad_params("Atmosphere Driver");
        auto& proc_params = ad_params.sublist("Atmosphere Processes");
        proc_params.set("Number of Entries", 1);
        proc_params.set<std::string>("Schedule Type", "Sequential");
        auto& p0 = proc_params.sublist("Process 0");
        p0.set<std::string>("Process Name", "RRTMGP");
        p0.set<std::string>("Grid","Physics");

        auto& gm_params = ad_params.sublist("Grids Manager");
        gm_params.set<std::string>("Type","User Provided");
        gm_params.set<std::string>("Reference Grid","Physics");

        // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
        // which rely on factory for process creation. The initialize method of the AD does that.
        auto& proc_factory = AtmosphereProcessFactory::instance();
        proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);

        // Need to register grids managers before we create the driver
        auto& gm_factory = GridsManagerFactory::instance();
        gm_factory.register_product("User Provided",create_user_provided_grids_manager);

        // Set the dummy grid in the UserProvidedGridManager
        // Recall that this class stores *static* members, so whatever
        // we set here, will be reflected in the GM built by the factory.
        UserProvidedGridsManager upgm;
        upgm.set_grid(std::make_shared<DummyPhysicsGrid>(num_cols));
        upgm.set_reference_grid("Physics");

        // Create a MPI communicator
        ekat::Comm atm_comm (MPI_COMM_WORLD);

        // Create the driver
        AtmosphereDriver ad;

        // Dummy timestamp
        util::TimeStamp time (0,0,0,0);

        // Initialize the driver, run the driver, cleanup
        ad.initialize(atm_comm, ad_params, time);
        ad.run(300.0);
        ad.finalize();
        upgm.clean_up();

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
