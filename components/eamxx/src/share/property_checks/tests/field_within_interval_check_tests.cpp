#include <catch2/catch.hpp>

#include "pc_tests_helpers.hpp"

#include "share/property_checks/field_within_interval_check.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("field_within_interval_check") {
  using namespace ShortFieldTagsNames;

  auto seed = get_random_test_seed();

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_lcols = 2;
  const int nlevs = 12;

  // Create grid
  auto grid = create_test_grid(comm,num_lcols,nlevs);

  // Create field to test, and extra data field
  auto f = create_test_field(grid);
  auto data = create_data_field(grid);

  // Check that the values of a field lie within an interval.
  SECTION ("without_repair_bounds") {

    // Create the check
    auto interval_check = std::make_shared<FieldWithinIntervalCheck>(f, grid, 0, 1, true);
    REQUIRE(interval_check->can_repair());

    interval_check->set_additional_data_field(data);

    // Assign in-bound values to the field and make sure it passes the check
    randomize_uniform(f,seed++,0.01,0.99);
    auto res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Assign out-of-bounds values to the field, make sure it fails the check
    auto f_012 = f.subfield(COL,0).subfield(CMP,1).subfield(LEV,2);
    auto f_123 = f.subfield(COL,1).subfield(CMP,2).subfield(LEV,3);
    f_012.deep_copy(-1);
    f_123.deep_copy(2);
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);
    std::string expected_msg =
      "Check failed.\n"
      "  - check name: field_1 within interval [0, 1]\n"
      "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n"
      "  - minimum:\n"
      "    - value: -1\n"
      "    - indices (w/ global column index): (0,1,2)\n"
      "    - lat/lon: (0, 0)\n"
      "  - maximum:\n"
      "    - value: 2\n"
      "    - indices (w/ global column index): (1,2,3)\n"
      "    - lat/lon: (1, -1)\n";

    REQUIRE(res_and_msg.msg == expected_msg);

    // Check that after repairing the field, the check passes
    interval_check->repair();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);
  }

  SECTION ("with_repair_bounds") {
    using FWIC = FieldWithinIntervalCheck;

    // Repairable bounds must be less tight than bounds
    REQUIRE_THROWS (FWIC(f,grid,0,1,true,1,2));
    REQUIRE_THROWS (FWIC(f,grid,0,1,true,-1,0));

    // Create check
    auto interval_check = std::make_shared<FWIC>(f, grid, 0, 1, true,-1,2);
    REQUIRE(interval_check->can_repair());

    // Assign slightly out-of-bound values to the field
    randomize_uniform(f,seed++,0.01,0.99);
    auto f_beg = f.subfield(COL,0).subfield(CMP,0).subfield(LEV,0);
    auto f_end = f.subfield(COL,num_lcols-1).subfield(CMP,2).subfield(LEV,nlevs-1);
    f_beg.deep_copy(-0.5);
    f_end.deep_copy(1.5);
    auto res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);

    // Assign VERY out-of-bounds values to the field, make sure it fails the check
    f_beg.deep_copy(-2);
    f_end.deep_copy(3);
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Fail);
    std::string expected_msg =
      "Check failed.\n"
      "  - check name: field_1 within interval [0, 1]\n"
      "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n"
      "  - minimum:\n"
      "    - value: -2\n"
      "    - indices (w/ global column index): (0,0,0)\n"
      "    - lat/lon: (0, 0)\n"
      "  - maximum:\n"
      "    - value: 3\n"
      "    - indices (w/ global column index): (1,2,11)\n"
      "    - lat/lon: (1, -1)\n";

    REQUIRE(res_and_msg.msg == expected_msg);

    // Verify that, upon repairing the field, the check now passes
    interval_check->repair();
    f.sync_to_host();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Re-assign an out-of-bounds value to the field, but only in the lower-bound.  Check if it
    // fails and reports just the min fail.
    std::vector<int> exp_fail_loc = {0,0,0};
    f_beg.deep_copy(-2);
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Fail);
    REQUIRE(res_and_msg.fail_loc_indices == exp_fail_loc);
    // Repair for next check.
    interval_check->repair();
    f.sync_to_host();

    // Re-assign an out-of-bounds value to the field, but only in the upper-bound.  Check if it
    // fails and reports just the max fail.
    exp_fail_loc = {num_lcols-1,2,nlevs-1};
    f_end.deep_copy(4);
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Fail);
    REQUIRE(res_and_msg.fail_loc_indices == exp_fail_loc);
  }
}

} // namespace scream
