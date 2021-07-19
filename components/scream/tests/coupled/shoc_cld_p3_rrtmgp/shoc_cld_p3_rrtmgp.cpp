#include <catch2/catch.hpp>

// Boiler plate, needed for all runs
#include "control/atmosphere_driver.hpp"
#include "share/atm_process/atmosphere_process.hpp"
// Boiler plate, needed for when physics is part of the run
#include "physics/register_physics.hpp"
#include "physics/share/physics_only_grids_manager.hpp"
// EKAT headers
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

TEST_CASE("shoc-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  constexpr int num_iters = 10;

  // Initialize yakl
  if(!yakl::isInitialized()) { yakl::init(); }

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  register_physics();

  // Create the grids manager
  auto& gm_params = ad_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run 
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  for (int i=0; i<num_iters; ++i) {
    printf("Run iteration %2d\n",i);
    ad.run(300.0);
  }

  // Finalize 
  ad.finalize();
  // Finalize YAKL
  yakl::finalize();

  // If we got here, we were able to run shoc
  REQUIRE(true);
}

} // empty namespace
