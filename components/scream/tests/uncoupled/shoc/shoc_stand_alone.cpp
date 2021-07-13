#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/shoc/atmosphere_macrophysics.hpp"

#include "physics/share/physics_only_grids_manager.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

// === A dummy physics grids for this test === //

TEST_CASE("shoc-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  constexpr int num_iters = 10;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("shoc",&create_atmosphere_process<SHOCMacrophysics>);
  gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);

  // Create the grids manager
  auto& gm_params = ad_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  // Set pref_mid for shoc.  Note this is only temporary.
  // TODO: have pref_mid be part of the grid information and
  // have the shoc interface grab it from the grid, not the FM.
  const auto& grid = ad.get_grids_manager()->get_grid("Physics");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());
  int nlay = grid->get_num_vertical_levels();
  auto d_pref_mid = field_mgr.get_field("pref_mid").get_reshaped_view<Real*>();
  auto h_pref_mid = Kokkos::create_mirror_view(d_pref_mid);
  Kokkos::deep_copy(h_pref_mid,d_pref_mid);
  for (int k=0;k<nlay;k++) {
    h_pref_mid(k) = 1e5 - k*(1e5-8e5)/(nlay-1);
  }
  Kokkos::deep_copy(d_pref_mid,h_pref_mid);

  for (int i=0; i<num_iters; ++i) {
    ad.run(300.0);
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of SHOC, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run shoc
  REQUIRE(true);
}

} // empty namespace
