#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/register_physics.hpp"
#include "physics/mam/eamxx_mam_aci_process_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/logging/ekat_logger.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("mam4-aci-standalone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  ekat::logger::Logger<> logger("mam4-aci",
                                ekat::logger::LogLevel::debug, atm_comm);


  // Load ad parameter list
  std::string fname = "input_aci.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  logger.debug("aci yaml parsed.");

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  logger.info("running MAMAci standalone test with dt = {} for {} steps.", dt, nsteps);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  register_physics();
  register_mesh_free_grids_manager();
  register_diagnostics();
  logger.debug("products registered.");

  // Create the driver
  AtmosphereDriver ad;
  logger.debug("driver created.");

  // Init and run
  ad.initialize(atm_comm,ad_params,t0);
  logger.debug("driver initialized.");

  logger.info("Start time stepping loop ... [0%]");
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    logger.info(" Iteration {} completed; [{}]", i+1, 100*(i+1)/nsteps);
  }

  // Finalize
  ad.finalize();

  // If we got here, we were able to run mam4 aci
  REQUIRE(true);
}
} // namespace scream
