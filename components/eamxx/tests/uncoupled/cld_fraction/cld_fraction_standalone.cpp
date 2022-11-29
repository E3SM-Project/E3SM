#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/cld_fraction/atmosphere_cld_fraction.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("cld_fraction-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

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

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  ad.initialize(atm_comm,ad_params,t0);

  // Because this is a relatively simple test based on two variables, we initialize them here
  // rather than use the netCDF input structure.
  const auto& grid = ad.get_grids_manager()->get_grid("Point Grid");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());
  
  int num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  int num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  const auto& qi_field           = field_mgr.get_field("qi");
  const auto& liq_cld_frac_field = field_mgr.get_field("cldfrac_liq");  //TODO: This FM name will probably change soon.
  const auto& qi           = qi_field.get_view<Real**,Host>();
  const auto& liq_cld_frac = liq_cld_frac_field.get_view<Real**,Host>();

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
  // Sync the fields back to device
  qi_field.sync_to_dev();
  liq_cld_frac_field.sync_to_dev();

  // Run the code
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

  // Check ice and total cloud fraction values
  // Sync the values on device back to the host view.
  const auto& ice_cld_frac_field = field_mgr.get_field("cldfrac_ice");  //TODO: This FM name will probably change soon.
  const auto& tot_cld_frac_field = field_mgr.get_field("cldfrac_tot");   //TODO: This FM name will probably change soon.
  ice_cld_frac_field.sync_to_host();
  tot_cld_frac_field.sync_to_host();
  const auto& ice_cld_frac = ice_cld_frac_field.get_view<Real**,Host>();
  const auto& tot_cld_frac = tot_cld_frac_field.get_view<Real**,Host>();
  
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
      Real y_cmp = (std::sin(xval-phase)+1.0)/2.0;
      REQUIRE(qi(icol,jlev)==y_cmp);
      y_cmp = (std::sin(xval+phase)+1.0)/2.0;
      REQUIRE(liq_cld_frac(icol,jlev)==y_cmp);
      // Test that the cloud fraction calculation appropriately calculated the ice cloud fraction
      if (qi(icol,jlev)>1e-5)  // TODO: This may also be updated to be a tunable parameter, we might need to grab this value from the parameter list.
      {
        REQUIRE(ice_cld_frac(icol,jlev)==1.0);
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

}

} // empty namespace
