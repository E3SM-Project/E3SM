#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/spa/atmosphere_prescribed_aerosol.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("spa-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

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

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("SPA",&create_atmosphere_process<SPA>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);
  register_diagnostics();

  // Create the grids manager
  auto& gm_params = ad_params.sublist("grids_manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

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
}

} // empty namespace
