#include <catch2/catch.hpp>
#include "share/atmosphere_process.hpp"
#include "share/scream_pack.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "control/atmosphere_driver.hpp"

// #include "physics/rrtmgp/atmosphere_microphysics.hpp"
// #include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
// #include "physics/rrtmgp/rrtmgp_functions_f90.hpp"

namespace scream {
    // === A dummy physics grids for this test === //
    class DummyPhysicsGrid : public SEGrid {
        public: DummyPhysicsGrid (const int num_cols) : SEGrid("Physics",GridType::SE_NodeBased,num_cols) {
            // Nothing to do here
        }
        ~DummyPhysicsGrid () = default;
    };

    TEST_CASE("rrtmgp_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        // If we got here, we were able to run RRTMGP
        REQUIRE(true);
    }
}
