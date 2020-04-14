#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "share/scream_parameter_list.hpp"
#include "dummy_atm_setup.hpp"

namespace scream {

TEST_CASE ("dag_check","[!throws]")
{
  // This test checks that unmet dependencies in the Atm DAG will throw
  ParameterList ad_params("Atmosphere Driver");
  auto& proc_params = ad_params.sublist("Atmosphere Processes");

  proc_params.set("Number of Entries",1);
  proc_params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = proc_params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Physics_fwd");
  p0.set<int>("Number of vector components",1);

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Reference Grid","Physics_fwd");
  gm_params.set<std::string>("Type","User Provided");

  // Setup the atm factories and grid manager
  dummy_atm_init(1);

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  control::AtmosphereDriver ad;

  util::TimeStamp init_time(0,0,0);

  // Since Physics_fwd has an unmet dependency, this should throw
  REQUIRE_THROWS(ad.initialize(atm_comm,ad_params,init_time));

  // Cleanup atm factories and grids manager
  dummy_atm_cleanup();
}

} // namespace scream
