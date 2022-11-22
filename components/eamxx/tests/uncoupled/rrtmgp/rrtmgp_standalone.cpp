#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/rrtmgp/atmosphere_radiation.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("rrtmgp-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string inputfile = ekat::TestSession::get().params.at("inputfile");
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(inputfile,ad_params) );

  // Time stepping parameters
  auto& ts = ad_params.sublist("Time Stepping");
  const auto dt = ts.get<int>("Time Step");
  const auto start_date = ts.get<std::vector<int>>("Start Date");
  const auto start_time  = ts.get<std::vector<int>>("Start Time");
  const auto nsteps     = ts.get<int>("Number of Steps");

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

  util::TimeStamp t0 (start_date, start_time);
  EKAT_ASSERT_MSG (t0.is_valid(), "Error! Invalid start date.\n");

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("rrtmgp",&create_atmosphere_process<RRTMGPRadiation>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);
  register_diagnostics();

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  ad.initialize(atm_comm,ad_params,t0);

  // Get a pointer to the field manager so we can query fields
  const auto& grid = ad.get_grids_manager()->get_grid("Point Grid");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());

  // Get field managed variables we need to check
  auto sw_flux_up = field_mgr.get_field("sw_flux_up");

  // Create deep copies so that we can check values before and after call to ad.run
  // Note: we have to do some trickery here to make sure the new fields we allocate
  // get the same size as the (packed) ones in the FM, otherwise the vertical dim
  // might be padded in the FM fields and unpadded in our copies, which will cause
  // the deep_copy below to fail.
  Field sw_flux_up_old = Field(sw_flux_up.get_header().get_identifier());
  const auto& ap_new = sw_flux_up.get_header().get_alloc_properties();
  auto& ap_old = sw_flux_up_old.get_header().get_alloc_properties();
  ap_old.request_allocation(ap_new.get_largest_pack_size());
  sw_flux_up_old.allocate_view();

  // Start stepping
  if (atm_comm.am_i_root()) {
    printf("Start time stepping loop...       [  0%%]\n");
  }
  for (int i=0; i<nsteps; ++i) {

    // Create a (deep) copy of fields we want to check before calling ad.run() so we can verify
    // that these fields do or do not change as we expect them to based on the rad frequency
    sw_flux_up_old.deep_copy(sw_flux_up);

    ad.run(dt);
    if (atm_comm.am_i_root()) {
      std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
      std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
    }

    // Test that in between rad steps, we maintain the same values of fluxes and heating rates
    // get rad fluxes and heating rates before; we set rad_requency to 3 in the input.yaml, so
    // the first two steps should look the same
    auto d_sw_flux_up_new = sw_flux_up.get_view<Real**,Host>();
    auto d_sw_flux_up_old = sw_flux_up_old.get_view<Real**,Host>();
    if (i == 0) {
        REQUIRE(!views_are_equal(sw_flux_up_old, sw_flux_up));
    } else if (i == 1) {
        REQUIRE(views_are_equal(sw_flux_up_old, sw_flux_up));
    } else if (i == 2) {
        REQUIRE(views_are_equal(sw_flux_up_old, sw_flux_up));
    } else if (i == 3) {
        REQUIRE(!views_are_equal(sw_flux_up_old, sw_flux_up));
    }

  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of SHOC, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run shoc
  REQUIRE(true);
}

} // empty namespace
