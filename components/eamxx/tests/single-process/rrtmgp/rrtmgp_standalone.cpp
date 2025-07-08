#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "physics/register_physics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
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
  parse_yaml_file(inputfile,ad_params);

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

  // Need to register products in the factory *before* we create any atm process or grids manager.
  register_physics();
  register_mesh_free_grids_manager();
  register_diagnostics();

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  ad.initialize(atm_comm,ad_params,t0);

  // Get a pointer to the field manager so we can query fields
  const auto& grid = ad.get_grids_manager()->get_grid("point_grid");
  const auto& field_mgr = *ad.get_field_mgr();

  // Get field managed variables we need to check
  auto rad_heating_pdel = field_mgr.get_field("rad_heating_pdel", grid->name());

  // Create deep copies so that we can check values before and after call to ad.run
  // Note: we have to do some trickery here to make sure the new fields we allocate
  // get the same size as the (packed) ones in the FM, otherwise the vertical dim
  // might be padded in the FM fields and unpadded in our copies, which will cause
  // the deep_copy below to fail.
  Field rad_heating_pdel_old = Field(rad_heating_pdel.get_header().get_identifier());
  const auto& ap_new = rad_heating_pdel.get_header().get_alloc_properties();
  auto& ap_old = rad_heating_pdel_old.get_header().get_alloc_properties();
  ap_old.request_allocation(ap_new.get_largest_pack_size());
  rad_heating_pdel_old.allocate_view();

  int rad_freq = ad_params.sublist("atmosphere_processes").sublist("rrtmgp").get<int>("rad_frequency");
  // Start stepping
  if (atm_comm.am_i_root()) {
    printf("Start time stepping loop...       [  0%%]\n");
  }

  auto time = t0;
  for (auto time=t0+dt; time.get_num_steps()<=nsteps; time+=dt) {
    auto istep = time.get_num_steps();
    // Create a (deep) copy of fields we want to check before calling ad.run() so we can verify
    // that these fields do or do not change as we expect them to based on the rad frequency
    rad_heating_pdel_old.deep_copy(rad_heating_pdel);

    ad.run(dt);
    if (atm_comm.am_i_root()) {
      std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << istep << " completed";
      std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(istep)/nsteps << "%]\n";
    }

    // Test that in between rad steps, we maintain the same value of heating rate
    // get rad heating rates before; we set rad_requency to 3 in the input.yaml, so
    // the first two steps should look the same
    auto d_rad_heating_pdel_new = rad_heating_pdel.get_view<Real**,Host>();
    auto d_rad_heating_pdel_old = rad_heating_pdel_old.get_view<Real**,Host>();

    if (istep==1 or istep%rad_freq==0) {
      REQUIRE(!views_are_equal(rad_heating_pdel_old, rad_heating_pdel));
    } else {
      REQUIRE(views_are_equal(rad_heating_pdel_old, rad_heating_pdel));
    }
  }

  // Finalize
  ad.finalize();
}

} // empty namespace
