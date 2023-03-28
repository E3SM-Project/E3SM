#include "catch2/catch.hpp"

// The AD
#include "control/atmosphere_driver.hpp"

// Physics includes
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"
#include "physics/cld_fraction/atmosphere_cld_fraction.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "diagnostics/register_diagnostics.hpp"

// Dynamics includes
#include "dynamics/register_dynamics.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

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
  parse_yaml_file(fname,ad_params);
  ad_params.print();

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  // Register all atm procs and the grids manager in the respective factories
  register_dynamics();
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
  proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);
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
