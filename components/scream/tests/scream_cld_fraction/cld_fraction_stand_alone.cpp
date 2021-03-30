#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/cld_fraction/atmosphere_cld_fraction.hpp"

#include "physics/share/physics_only_grids_manager.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

TEST_CASE("cld_fraction-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  constexpr int num_iters = 1;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.,
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
  gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);

  // Create the grids manager
  auto& gm_params = ad_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run (do not finalize, or you'll clear the field repo!)
  util::TimeStamp time (0,0,0,0);

  ad.initialize(atm_comm,ad_params,time);

  // Because this is a relatively simple test based on two variables, we initialize them here
  // rather than use the netCDF input structure.
  auto& field_repo = ad.get_field_repo();
  const auto& grid = ad.get_grids_manager()->get_grid("Physics");
  
  int num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  int num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  const auto& qi = field_repo.get_field("qi","Physics").get_reshaped_view<Real**,Host>();
  const auto& liq_cld_frac = field_repo.get_field("alst","Physics").get_reshaped_view<Real**,Host>();  //TODO: This FM name will probably change soon.

  for (int icol=0;icol<num_cols;++icol)
  {
    for (int jlev=0;jlev<num_levs;++jlev)
    {
      Real phase = icol*3.14/2.0/num_cols;
      Real xval  = jlev*3.14/2.0/num_levs;
      // Assign a simple functional value to qi based on the column index and level
      qi(icol,jlev) = (1.0+std::sin(xval-phase))/2.0;
      // Do something similar for liq_cloud_fraction
      liq_cld_frac(icol,jlev) = (1.0+std::sin(xval+phase))/2.0;
    }
  }

  // Run the code
  for (int i=0; i<num_iters; ++i) {
    ad.run(300.0);
  }

  // Check ice and total cloud fraction values
  // Sync the values on device back to the host view.
  const auto& ice_cld_frac_field = field_repo.get_field("aist","Physics");
  const auto& tot_cld_frac_field = field_repo.get_field("ast","Physics");
  const auto& ice_cld_frac_dev = ice_cld_frac_field.get_reshaped_view<Real**>();  //TODO: This FM name will probably change soon.
  const auto& tot_cld_frac_dev = tot_cld_frac_field.get_reshaped_view<Real**>();  //TODO: This FM name will probably change soon.
  auto ice_cld_frac = Kokkos::create_mirror_view(ice_cld_frac_dev);
  auto tot_cld_frac = Kokkos::create_mirror_view(tot_cld_frac_dev);
  Kokkos::deep_copy(ice_cld_frac,ice_cld_frac_dev);
  Kokkos::deep_copy(tot_cld_frac,tot_cld_frac_dev);
  
  for (int icol=0;icol<num_cols;++icol)
  {
    for (int jlev=0;jlev<num_levs;++jlev)
    {
      // Make sure there are no cloud fractions greater than 1.0 or less than 0.0
      REQUIRE((!(liq_cld_frac(icol,jlev)>1.0) or !(ice_cld_frac(icol,jlev)>1.0) or !(tot_cld_frac(icol,jlev)>1.0)));
      REQUIRE((!(liq_cld_frac(icol,jlev)<0.0) or !(ice_cld_frac(icol,jlev)<0.0) or !(tot_cld_frac(icol,jlev)<0.0)));
      // make sure that the cloud fraction calculation didn't accidentally change qi or liq_cld_fraction
      Real phase = icol*3.14/2.0/num_cols;
      Real xval  = jlev*3.14/2.0/num_levs;
      REQUIRE(qi(icol,jlev)==(std::sin(xval-phase)+1.0)/2.0);
      REQUIRE(liq_cld_frac(icol,jlev)==(std::sin(xval+phase)+1.0)/2.0);
      // Test that the cloud fraction calculation appropriately calculated the ice cloud fraction
      if (qi(icol,jlev)>1e-5)  // TODO: This may also be updated to be a tunable parameter, we might need to grab this value from the parameter list.
      {
        //REQUIRE(ice_cld_frac(icol,jlev)==1.0);
        EKAT_REQUIRE_MSG(ice_cld_frac(icol,jlev)==1.0, "error - " + std::to_string(icol) + " " + std::to_string(jlev) + " " + std::to_string(qi(icol,jlev)));
      } else {
        REQUIRE(ice_cld_frac(icol,jlev)==0.0);
      }
      // Test that the total cloud fraction is correctly calculated
      if (liq_cld_frac(icol,jlev) >= ice_cld_frac(icol,jlev))
      {
        REQUIRE(tot_cld_frac(icol,jlev)==liq_cld_frac(icol,jlev));
      } else {
        REQUIRE(tot_cld_frac(icol,jlev)==ice_cld_frac(icol,jlev));
      }
    }
  }

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run cld_fraction
  REQUIRE(true);
}

} // empty namespace
