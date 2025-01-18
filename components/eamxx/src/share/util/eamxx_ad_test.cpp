#include "catch2/catch.hpp"

// The AD
#include "control/atmosphere_driver.hpp"

// Physcis/dynamics/diagnostic includes
#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "dynamics/register_dynamics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

// EKAT headers
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

/*
 * The atm configuration created in this simple test is fully
 * configurable from input yaml file, just like any atm instance.
 * Notice that, despite the fact that the source code is the same
 * for all atm configurations, in order to use a particular atm
 * proc (e.g., a physics package), this source code must be linked
 * against the lib providing such atm proc, otherwise it won't
 * get registered in the atm proc factory.
 */

TEST_CASE("scream_ad_test") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // User can prescribe input file name via --args --ifile <file>
  auto& session = ekat::TestSession::get();
  session.params.emplace("ifile","input.yaml");
  std::string fname = session.params["ifile"];

  // Load ad parameter list
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  ad_params.print();

  // Time stepping parameters
        auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  run_t0_str = ts.get<std::string>("run_t0");
  const auto  run_t0     = util::str_to_time_stamp(run_t0_str);
  const auto  case_t0_str = ts.get<std::string>("case_t0",run_t0_str);
  const auto  case_t0     = util::str_to_time_stamp(case_t0_str);

  // Register all atm procs, grids manager, and diagnostics in the respective factories
  register_dynamics();
  register_physics();
  register_diagnostics();
  register_mesh_free_grids_manager();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
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
}
