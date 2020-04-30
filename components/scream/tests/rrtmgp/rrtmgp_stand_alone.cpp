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

// #include "physics/rrtmgp/atmosphere_microphysics.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
// #include "physics/rrtmgp/rrtmgp_functions_f90.hpp"

namespace scream {
    // === A dummy physics grids for this test === //
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

        // Do something interesting here...
        // NOTE: these will get replaced with AD stuff that handles these
        rrtmgp_init_f90();
        rrtmgp_main_f90();
        rrtmgp_finalize_f90();

        // If we got here, we were able to run the above code
        REQUIRE(true);
    }
}
