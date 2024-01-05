#include "catch2/catch.hpp"

// The driver
#include "control/atmosphere_driver.hpp"

// DYNAMICS and PHYSICS includes
#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "physics/shoc/eamxx_shoc_process_interface.hpp"
#include "physics/cld_fraction/eamxx_cld_fraction_process_interface.hpp"
#include "physics/rrtmgp/eamxx_rrtmgp_process_interface.hpp"
#include "dynamics/register_dynamics.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "diagnostics/register_diagnostics.hpp"

// EKAT headers
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_assert.hpp"

// Hommexx includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"
#include "FunctorsBuffersManager.hpp"

TEST_CASE("scream_homme_physics", "scream_homme_physics") {
  using namespace scream;
  using namespace scream::control;

  // Load ad parameter list
  const auto& session = ekat::TestSession::get();
  std::string fname = session.params.at("ifile");
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  // NOTE: we register physics by hand (rather than with register_physics()), since
  //       we don't want to require spa to run this test
  register_dynamics();
  register_diagnostics();

  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
  proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  run_t0_str  = ts.get<std::string>("run_t0");
  const auto  run_t0      = util::str_to_time_stamp(run_t0_str);
  const auto  case_t0_str = ts.get<std::string>("case_t0");
  const auto  case_t0     = util::str_to_time_stamp(case_t0_str);

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  // NOTE: Kokkos is finalize in ekat_catch_main.cpp, and YAKL is finalized
  //       during RRTMGPRatiation::finalize_impl, after RRTMGP has deallocated
  //       all its arrays.
  ad.initialize(atm_comm,ad_params,run_t0,case_t0);

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


  // If we got here, we were able to run homme
  REQUIRE(true);
}
