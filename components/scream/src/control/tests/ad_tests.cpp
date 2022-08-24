#include "dummy_atm_setup.hpp"

#include "control/atmosphere_driver.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <catch2/catch.hpp>

namespace scream {

TEST_CASE ("ad_tests","[!throws]")
{
  // Load ad parameter list
  std::string fname = "ad_tests.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Setup the atm factories and grid manager
  dummy_atm_init();

  // Create the driver
  control::AtmosphereDriver ad;

  // Init and run a single time step
  util::TimeStamp t0(2000,1,1,0,0,0);
  ad.initialize(atm_comm,ad_params,t0);

  // Verify that the atm proc group has the expected specs
  const auto& apg = ad.get_atm_processes();
  REQUIRE (apg->get_num_processes()==3);
  REQUIRE (apg->get_fields_in().size()==3);
  REQUIRE (apg->get_groups_in().size()==0);
  REQUIRE (apg->get_fields_out().size()==3);
  REQUIRE (apg->get_groups_out().size()==1);

  // Check correct initialization of the input fields for the 1st process
  const auto& f_in = apg->get_fields_in();
  for (const auto& f : f_in) {
    const auto& fn = f.get_header().get_identifier().name();

    // Create 'target' field
    Field check(f.get_header().get_identifier());
    check.allocate_view();

    // Fill target field based on what IC were in the yaml file
    if (fn=="A" || fn=="Z") {
      check.deep_copy(1.0);
    } else if (fn=="V") {
      check.get_component(0).deep_copy(2.0);
      check.get_component(1).deep_copy(3.0);
    } else {
      EKAT_ERROR_MSG ("Error! Unexpected input field for this test.\n");
    }

    // Check the field
    REQUIRE (views_are_equal(f,check));
  }

  // Run ad
  ad.run(10);

  // At this point, output fields should have timestamp updated
  const auto& f_out = apg->get_fields_out();
  auto t = t0 + 10;
  for (const auto& f : f_out) {
    const auto& f_ts = f.get_header().get_tracking().get_time_stamp();
    REQUIRE (f_ts==t);
  }

  // Cleanup
  ad.finalize ();
  dummy_atm_cleanup();
}

} // namespace scream
