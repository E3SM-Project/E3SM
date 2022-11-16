#include "catch2/catch.hpp"

// The AD
#include "control/atmosphere_driver.hpp"

// Dynamics includes
#include "dynamics/register_dynamics.hpp"

// Physics includes
#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"

// EKAT headers
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_assert.hpp"

TEST_CASE("scream_homme_physics", "scream_homme_physics") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Time stepping parameters
  ad_params.print();
  auto& ts = ad_params.sublist("Time Stepping");
  const auto dt = ts.get<int>("Time Step");
  const auto start_date = ts.get<std::vector<int>>("Start Date");
  const auto start_time = ts.get<std::vector<int>>("Start Time");
  const auto nsteps     = ts.get<int>("Number of Steps");

  util::TimeStamp t0 (start_date, start_time);
  EKAT_ASSERT_MSG (t0.is_valid(), "Error! Invalid start date.\n");

  // Register all atm procs and the grids manager in the respective factories
  register_dynamics();
  register_physics();
  register_diagnostics();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  // NOTE: Kokkos is finalize in ekat_catch_main.cpp, and YAKL is finalized
  //       during RRTMGPRatiation::finalize_impl, after RRTMGP has deallocated
  //       all its arrays.
  ad.initialize(atm_comm,ad_params,t0);

  if (atm_comm.am_i_root()) {
    printf("Start time stepping loop...       [  0%%]\n");
  }
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    if (atm_comm.am_i_root()) {
      std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
      std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
    }
  }
  ad.finalize();


  // If we got here, we were able to run without errors.
  REQUIRE (true);
}
