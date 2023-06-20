#include "catch2/catch.hpp"

// The AD
#include "control/atmosphere_driver.hpp"

// Dynamics includes
#include "dynamics/register_dynamics.hpp"

// Physics includes
#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "physics/mam/eamxx_mam_microphysics_process_interface.hpp"

// EKAT headers
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_assert.hpp"

TEST_CASE("scream_homme_physics", "scream_homme_physics_mam4") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  ekat::logger::Logger<> logger("homme-mam4",
                                ekat::logger::LogLevel::debug, atm_comm);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  ad_params.print();

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  logger.info("running HOMME/MAMMicrophysics coupled test with dt = {} for {} steps.", dt, nsteps);

  // Register all atm procs and the grids manager in the respective factories
  register_dynamics();
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("MAMMicrophysics",&create_atmosphere_process<MAMMicrophysics>);
//   logger.debug("registered products: {}", proc_factory.print_registered_products());
// TODO: register_diagnostics();


  // Create the driver
  AtmosphereDriver ad;
  logger.debug("driver created.");

  // Init, run, and finalize
  // NOTE: Kokkos is finalize in ekat_catch_main.cpp, and YAKL is finalized
  //       during RRTMGPRatiation::finalize_impl, after RRTMGP has deallocated
  //       all its arrays.
  ad.initialize(atm_comm,ad_params,t0);
  logger.debug("driver initialized.");

  logger.info("Start time stepping loop ... [0%]");
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    logger.info(" Iteration {} completed; [{}]", i+1, 100*(i+1)/nsteps);
  }
  ad.finalize();


  // If we got here, we were able to run without errors.
  REQUIRE (true);
}
