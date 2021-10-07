#include "dummy_atm_setup.hpp"

#include "control/atmosphere_driver.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <catch2/catch.hpp>

namespace scream {

TEST_CASE ("group_requirements","[!throws]")
{
  constexpr int num_cols = 4;
  constexpr int num_vl   = 2;

  // Load ad parameter list
  std::string fname = "ad_tests.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Setup the atm factories and grid manager
  dummy_atm_init(num_cols, num_vl, atm_comm);

  // Create the driver
  control::AtmosphereDriver ad;

  // Init and run a single time step
  util::TimeStamp init_time(2000,1,1,0,0,0);
  ad.initialize(atm_comm,ad_params,init_time);

  // Verify that the atm proc group has 1 required field and 0 required groups.
  const auto& apg = ad.get_atm_processes();
  REQUIRE (apg->get_num_processes()==3);
  REQUIRE (apg->get_fields_in().size()==1);
  REQUIRE (apg->get_groups_in().size()==0);

  auto field_mgr = ad.get_ref_grid_field_mgr();
  // Resume ad run and testing.
  ad.run(10);
  ad.finalize ();

  // Cleanup atm factories and grids manager
  dummy_atm_cleanup();
}

} // namespace scream
