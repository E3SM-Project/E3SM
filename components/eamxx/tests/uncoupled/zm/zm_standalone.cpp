#include <catch2/catch.hpp>
#include "control/atmosphere_driver.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/zm/atmosphere_deep_convection.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

#include <iostream>
#include <iomanip>

namespace scream {

TEST_CASE("zm-standalone", "") {
  using namespace scream;
  using namespace scream::control;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  
  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  
  // Need to register grids managers before we create the driver
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("ZM",&create_atmosphere_process<ZMDeepConvection>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);
  register_diagnostics();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  ad.initialize(atm_comm,ad_params,t0);
  printf("Start time stepping loop...       [  0%%]\n");
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
    std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
  }
  ad.finalize();

  // If we got here, we were able to run zm 
  REQUIRE(true);
}

} // namespace scream
