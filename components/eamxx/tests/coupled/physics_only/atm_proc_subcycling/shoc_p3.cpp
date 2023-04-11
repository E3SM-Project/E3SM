#include <catch2/catch.hpp>

// Boiler plate, needed for all runs
#include "control/atmosphere_driver.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

// Physics headers
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

// EKAT headers
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("shoc-p3", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  const auto& session = ekat::TestSession::get();
  std::string fname = session.params.at("ifile");
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Time stepping parameters
  auto& ts = ad_params.sublist("Time Stepping");
  const auto dt = ts.get<int>("Time Step");
  const auto start_date = ts.get<std::vector<int>>("Start Date");
  const auto start_time = ts.get<std::vector<int>>("Start Time");
  const auto nsteps     = ts.get<int>("Number of Steps");

  util::TimeStamp t0 (start_date, start_time);
  EKAT_ASSERT_MSG (t0.is_valid(), "Error! Invalid start date.\n");

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);
  register_mesh_free_grids_manager();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
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

  // If we got here, we were able to run shoc
  REQUIRE(true);
}

} // empty namespace
