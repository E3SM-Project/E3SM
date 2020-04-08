#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"

// Hommexx includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"

// P3 includes
#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/scream_p3_interface.hpp"

TEST_CASE("scream_homme_p3", "scream_homme_p3") {
  using namespace scream;
  using namespace scream::control;

  std::cout << "Start Test...\n";

  // Create a parameter list for inputs
  std::cout << "Create parameter list...";
  ParameterList ad_params("Atmosphere Driver");
  auto& params = ad_params.sublist("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "dynamics");
  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "P3");
  p1.set<std::string>("Grid","SE Physics");

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Type","Dynamics Driven");
  gm_params.set<std::string>("Reference Grid","SE Physics");
  std::cout << "Done\n";

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  std::cout << "Create factory...";
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("dynamics",&create_atmosphere_process<HommeDynamics>);
  proc_factory.register_product("P3",&create_atmosphere_process<P3Microphysics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Dynamics Driven",create_dynamics_driven_grids_manager);
  std::cout << "Done\n";

  // Create a comm
  std::cout << "Create comm...";
  Comm atm_comm (MPI_COMM_WORLD);
  std::cout << "Done\n";

  // Create the driver
  std::cout << "Create driver...";
  AtmosphereDriver ad;
  std::cout << "Done\n";

  // Init, run, and finalize
  util::TimeStamp time (0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  // Have to wait till now to get homme's parameters, cause it only gets init-ed during the driver initialization
  const auto& sp = Homme::Context::singleton().get_simulation_params();
  const int nmax = get_homme_param_value<int>("nmax");
  const int num_dyn_iters = 10; //nmax / (sp.qsplit*sp.rsplit);
  const double dt = get_homme_param_value<Real>("dt");

  for (int i=0; i<num_dyn_iters; ++i) {
    ad.run(dt);
  }
  ad.finalize();

  // If we got here, we were able to run homme
  REQUIRE(true);
}
