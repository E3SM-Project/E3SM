#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "dynamics/register_dynamics.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "dynamics/homme/homme_dimensions.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_assert.hpp"

// Hommexx includes
#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "ElementsGeometry.hpp"
#include "TimeLevel.hpp"

#include <iomanip>

TEST_CASE("scream_homme_standalone", "scream_homme_standalone") {
  using namespace scream;
  using namespace scream::control;

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
  register_dynamics();
  register_diagnostics();

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  ad.initialize(atm_comm,ad_params,t0);

  // Check that topography data from the FM matches Homme.
  {
    auto& geo = Homme::Context::singleton().get<Homme::ElementsGeometry>();
    auto phis = geo.m_phis;

    const auto& atm_process_group = ad.get_atm_processes();
    const auto& process = atm_process_group->get_process(0);
    auto homme_process = std::dynamic_pointer_cast<const HommeDynamics>(process);
    EKAT_REQUIRE_MSG (process, "Error! Cast to HommeDynamics failed.\n");

    const auto phinh_i = homme_process->get_internal_field("phi_int_dyn","Dynamics").get_view<Real****>();

    int nelem = Homme::Context::singleton().get<Homme::ElementsGeometry>().num_elems();
    constexpr int NVL = HOMMEXX_NUM_PHYSICAL_LEV;

    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*NP*NP),
                         KOKKOS_LAMBDA (const int idx) {
      const int ie = idx/(NP*NP);
      const int ip = (idx/NP)%NP;
      const int jp = idx%NP;
      EKAT_KERNEL_ASSERT(phinh_i(ie,ip,jp,NVL) == phis(ie,ip,jp));
    });
  }

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
