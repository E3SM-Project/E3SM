#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

// P3 includes
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/scream_p3_interface.hpp"

TEST_CASE("scream_homme_p3", "scream_homme_p3") {
  using namespace scream;
  using namespace scream::control;

  std::cout << "Start Test...\n";

  // Load ad parameter list
  std::cout << "Loading parameter list ... ";
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );
  std::cout << "done!\n";

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  std::cout << "Create factory ... ";
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("dynamics",&create_atmosphere_process<HommeDynamics>);
  proc_factory.register_product("P3",&create_atmosphere_process<P3Microphysics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Dynamics Driven",create_dynamics_driven_grids_manager);
  std::cout << "done!\n";

  // Create a comm
  std::cout << "Create atm comm ... ";
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  std::cout << "done!\n";

  // Create the driver
  std::cout << "Create and init atm driver ... ";
  AtmosphereDriver ad;

  // Init, run, and finalize
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);
  std::cout << "done!\n";

  // Have to wait till now to get homme's parameters, cause it only gets init-ed during the driver initialization
  const int num_dyn_iters = 10; //nmax / (sp.qsplit*sp.rsplit);
  const double dt = get_homme_param<Real>("dt");

  std::cout << "Beginning time loop ...\n";
  for (int i=0; i<num_dyn_iters; ++i) {
    ad.run(dt);
  }

  std::cout << "Time loop complete.\n"
            << "Finalizing atm driver ...";
  ad.finalize();
  std::cout << "done!\n";

  // If we got here, we were able to run homme
  REQUIRE(true);
}
