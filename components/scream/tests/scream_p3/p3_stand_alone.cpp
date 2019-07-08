#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "control/surface_coupling.hpp"
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/scream_p3_interface.hpp"

TEST_CASE("scream_p3_stand_alone", "scream_p3_stand_alone") {
  using namespace scream;
  using namespace scream::control;

  // Create a parameter list for inputs
  ParameterList ad_params("Atmosphere Driver");
  auto& params = ad_params.sublist("Atmosphere Processes");

  params.set("Number of Entries",1);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "microphysics");

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Type","Physics Driven");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("microphysics",&create_atmosphere_process<P3Microphysics>);

//  // Need to register grids managers before we create the driver
//  auto& gm_factory = GridsManagerFactory::instance();
//  gm_factory.register_product("Physics Driven",create_physics_driven_grids_manager);

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  ad.initialize(atm_comm,ad_params);

  // Have to wait till now to get homme's parameters, cause it only gets init-ed during the driver initialization
  const int num_p3_iters = 10;

  for (int i=0; i<num_p3_iters; ++i) {
    ad.run();
  }
  ad.finalize();

  // If we got here, we were able to run p3
  REQUIRE(true);
}
