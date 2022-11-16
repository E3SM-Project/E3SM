#include "catch2/catch.hpp"

// The driver
#include "control/atmosphere_driver.hpp"

// DYNAMICS and PHYSICS includes
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"
#include "physics/cld_fraction/atmosphere_cld_fraction.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
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
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

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
  auto& ts_pl = ad_params.sublist("Time Stepping");
  const auto dt = ts_pl.get<int>("Time Step");
  const auto run_start_date = ts_pl.get<std::vector<int>>("Run Start Date");
  const auto run_start_time = ts_pl.get<std::vector<int>>("Run Start Time");
  const auto case_start_date = ts_pl.get<std::vector<int>>("Case Start Date");
  const auto case_start_time = ts_pl.get<std::vector<int>>("Case Start Time");

  util::TimeStamp run_t0 (run_start_date, run_start_time);
  util::TimeStamp case_t0 (case_start_date, case_start_time);

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
  const auto nsteps = ts_pl.get<int>("Number of Steps");
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
