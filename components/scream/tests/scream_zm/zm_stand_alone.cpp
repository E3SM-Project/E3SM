#include <catch2/catch.hpp>
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "control/atmosphere_driver.hpp"
#include "physics/share/physics_only_grids_manager.hpp"

#include "physics/zm/atmosphere_deep_convection.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <iostream>
namespace scream {

TEST_CASE("zm-standalone", "") {
  
  using namespace scream;
  using namespace scream::control;
  constexpr int num_iters = 20;
//  constexpr int num_cols  = 32;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  
  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  
  // Need to register grids managers before we create the driver
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("ZM",&create_atmosphere_process<ZMDeepConvection>);
  gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);


  // Create the grids manager
  auto& gm_params = ad_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);


  // Create the driver
  AtmosphereDriver ad;

  // Init and run (do not finalize, or you'll clear the field repo!)
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);
  for (int i=0; i<num_iters; ++i) {
    ad.run(300.0);
}

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run zm 
  REQUIRE(true);
}

} // empty namespace
