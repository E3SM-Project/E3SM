#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "dynamics/register_dynamics.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_feutils.hpp"
#include "ekat/ekat_assert.hpp"

// Hommexx includes
#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"

#include <iomanip>

static int get_default_fpes () {
#ifdef SCREAM_FPE
  return (FE_DIVBYZERO |
          FE_INVALID   |
          FE_OVERFLOW);
#else
  return 0;
#endif
}

TEST_CASE("scream_homme_standalone", "scream_homme_standalone") {
  using namespace scream;
  using namespace scream::control;

  ekat::enable_fpes(get_default_fpes());

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Time stepping parameters
  auto& ts = ad_params.sublist("Time Stepping");
  const auto dt = ts.get<int>("Time Step");
  const auto start_date = ts.get<std::vector<int>>("Start Date");
  const auto start_time  = ts.get<std::vector<int>>("Start Time");
  const auto nsteps     = ts.get<int>("Number of Steps");

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

  util::TimeStamp t0 (start_date, start_time);
  EKAT_ASSERT_MSG (t0.is_valid(), "Error! Invalid start date.\n");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("dynamics",&create_atmosphere_process<HommeDynamics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Dynamics Driven",create_dynamics_driven_grids_manager);

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  ad.initialize(atm_comm,ad_params,t0);

  // Add checks to verify AD memory buffer and Homme FunctorsBuffersManager
  // are the same size and reference the same memory.
  auto& fbm  = Homme::Context::singleton().get<Homme::FunctorsBuffersManager>();
  auto& memory_buffer = ad.get_memory_buffer();
  EKAT_ASSERT_MSG(fbm.allocated_size()*sizeof(Real) == (long unsigned int)memory_buffer.allocated_bytes(),
                  "Error! AD memory buffer and Homme FunctorsBuffersManager have mismatched sizes.");
  EKAT_ASSERT_MSG(fbm.get_memory() == memory_buffer.get_memory(),
                  "Error! AD memory buffer and Homme FunctorsBuffersManager reference different memory.");

  printf("Start time stepping loop...       [  0%%]\n");
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
    std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
  }
  ad.finalize();

  // If we got here, we were able to run homme
  REQUIRE(true);
}
