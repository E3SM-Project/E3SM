#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "ekat/scream_parameter_list.hpp"
#include "ekat/scream_parse_yaml_file.hpp"
#include "dummy_atm_setup.hpp"

namespace scream {

TEST_CASE ("dag_check","[!throws]")
{
  // This test checks that unmet dependencies in the Atm DAG will throw

  // Load ad parameter list
  std::string fname = "ad_tests.yaml";
  ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Setup the atm factories and grid manager
  dummy_atm_init(1);

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  control::AtmosphereDriver ad;

  util::TimeStamp init_time(0,0,0,0.0);

  // Since Physics_fwd has an unmet dependency, this should throw
  REQUIRE_THROWS(ad.initialize(atm_comm,ad_params,init_time));

  // Cleanup atm factories and grids manager
  dummy_atm_cleanup();
}

} // namespace scream
