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
  parse_yaml_file(fname,ad_params);

  // Time stepping parameters
  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);

  // Cloud fraction specific parameters
  const auto procs_params   = ad_params.sublist("atmosphere_processes");
  const auto cld_params     = procs_params.sublist("CldFraction");
  const auto ice_thresh     = cld_params.get<double>("ice_cloud_threshold");
  const auto ice_thresh_out = cld_params.get<double>("ice_cloud_for_analysis_threshold");

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

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
  const auto& liq_cld_frac_field = field_mgr.get_field("cldfrac_liq");
  const auto& qi                 = qi_field.get_view<Real**,Host>();
  const auto& liq_cld_frac       = liq_cld_frac_field.get_view<Real**,Host>();

  // Set initial ice mass mixing ratio and liquid cloud fraction
  const Real Pi = 3.14159265358979323846;
        Real qi_amp;
        Real cf_amp;
  for (int icol=0;icol<num_cols;++icol)
  {
    // ensure edge case where no `qi` is present in profile
    qi_amp = icol == num_cols-1 ? 0.0 : 1e-3;
    for (int jlev=0;jlev<num_levs;++jlev)
    {
      cf_amp = jlev == num_levs-1 ? 0.0 : 0.6;
      Real phase = icol*Pi/2.0/(num_cols-1);
      Real xval  = jlev*Pi/2.0/(num_levs-1);
      // Assign a simple functional value to qi based on the column index and level
      qi(icol,jlev) = qi_amp*(1.0-std::sin(xval-phase))/2.0;
      // Do something similar for liq_cloud_fraction
      liq_cld_frac(icol,jlev) = cf_amp*(1.0+std::sin(xval+phase))/2.0;
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
  const auto& ice_cld_frac_field = field_mgr.get_field("cldfrac_ice");
  const auto& tot_cld_frac_field = field_mgr.get_field("cldfrac_tot");
  ice_cld_frac_field.sync_to_host();
  tot_cld_frac_field.sync_to_host();
  const auto& ice_cld_frac = ice_cld_frac_field.get_view<Real**,Host>();
  const auto& tot_cld_frac = tot_cld_frac_field.get_view<Real**,Host>();

  const auto& ice_cld_frac_field_4out = field_mgr.get_field("cldfrac_ice_for_analysis"); 
  const auto& tot_cld_frac_field_4out = field_mgr.get_field("cldfrac_tot_for_analysis"); 
  ice_cld_frac_field_4out.sync_to_host();
  tot_cld_frac_field_4out.sync_to_host();
  const auto& ice_cld_frac_4out = ice_cld_frac_field_4out.get_view<Real**,Host>();
  const auto& tot_cld_frac_4out = tot_cld_frac_field_4out.get_view<Real**,Host>();
  
  for (int icol=0;icol<num_cols;++icol)
  {
    qi_amp = icol == num_cols-1 ? 0.0 : 1e-3;
    for (int jlev=0;jlev<num_levs;++jlev)
    {
      cf_amp = jlev == num_levs-1 ? 0.0 : 0.6;
      // Make sure there are no cloud fractions greater than 1.0 or less than 0.0
      REQUIRE((!(liq_cld_frac(icol,jlev)>1.0) or !(ice_cld_frac(icol,jlev)>1.0) or !(tot_cld_frac(icol,jlev)>1.0)));
      REQUIRE((!(liq_cld_frac(icol,jlev)<0.0) or !(ice_cld_frac(icol,jlev)<0.0) or !(tot_cld_frac(icol,jlev)<0.0)));
      // make sure that the cloud fraction calculation didn't accidentally change qi or liq_cld_fraction
      Real phase = icol*Pi/2.0/(num_cols-1);
      Real xval  = jlev*Pi/2.0/(num_levs-1);
      Real y_cmp = qi_amp*(1.0-std::sin(xval-phase))/2.0;
      REQUIRE(qi(icol,jlev)==y_cmp);
      y_cmp = cf_amp*(std::sin(xval+phase)+1.0)/2.0;
      REQUIRE(liq_cld_frac(icol,jlev)==y_cmp);
      // Test that the cloud fraction calculation appropriately calculated the ice cloud fraction
      if (qi(icol,jlev)>ice_thresh) 
      {
        REQUIRE(ice_cld_frac(icol,jlev)==1.0);
      } else {
        REQUIRE(ice_cld_frac(icol,jlev)==0.0);
      }
      if (qi(icol,jlev)>ice_thresh_out) 
      {
        REQUIRE(ice_cld_frac_4out(icol,jlev)==1.0);
      } else {
        REQUIRE(ice_cld_frac_4out(icol,jlev)==0.0);
      }
      // Test that the total cloud fraction is correctly calculated
      if (liq_cld_frac(icol,jlev) >= ice_cld_frac(icol,jlev))
      {
        REQUIRE(tot_cld_frac(icol,jlev)==liq_cld_frac(icol,jlev));
      } else {
        REQUIRE(tot_cld_frac(icol,jlev)==ice_cld_frac(icol,jlev));
      }
      if (liq_cld_frac(icol,jlev) >= ice_cld_frac_4out(icol,jlev))
      {
        REQUIRE(tot_cld_frac_4out(icol,jlev)==liq_cld_frac(icol,jlev));
      } else {
        REQUIRE(tot_cld_frac_4out(icol,jlev)==ice_cld_frac_4out(icol,jlev));
      }
    }
  }

  // Finalize 
  ad.finalize();

}

} // empty namespace
