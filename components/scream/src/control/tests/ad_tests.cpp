#include "dummy_atm_setup.hpp"

#include "control/atmosphere_driver.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <catch2/catch.hpp>

namespace scream {

TEST_CASE ("group_requirements","[!throws]")
{
  constexpr int num_cols = 4;
  constexpr int num_vl   = 2;

  // Load ad parameter list
  std::string fname = "ad_tests.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Setup the atm factories and grid manager
  dummy_atm_init(num_cols, num_vl, atm_comm);

  // Create the driver
  control::AtmosphereDriver ad;

  // Init and run a single time step
  util::TimeStamp init_time(0,0,0,0.0);
  ad.initialize(atm_comm,ad_params,init_time);
  ad.run(1.0);
  ad.finalize ();

  // Cleanup atm factories and grids manager
  dummy_atm_cleanup();
}

} // namespace scream
