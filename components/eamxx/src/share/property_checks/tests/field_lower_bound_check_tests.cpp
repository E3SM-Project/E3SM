#include <catch2/catch.hpp>

#include "pc_tests_helpers.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("lower_bound") {
  auto seed = get_random_test_seed();

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_lcols = 2;
  const int nlevs = 12;

  // Create grid
  auto grid = create_test_grid(comm,num_lcols,nlevs);

  // Create field to test, and extra data field
  auto f = create_test_field(grid);
  auto data = create_data_field(grid);
  auto ones = f.clone();
  ones.deep_copy(1);

  // Create the check
  auto lower_bound_check = std::make_shared<FieldLowerBoundCheck>(f,grid,-1.0, true);
  REQUIRE(lower_bound_check->can_repair());

  // Assign in-bound values to the field and make sure it passes the check
  randomize_uniform(f,seed,0.01,0.99);
  auto res_and_msg = lower_bound_check->check();
  REQUIRE(res_and_msg.result==CheckResult::Pass);

  // Assign out-of-bounds values to the field, make sure it fails the check
  f.update(ones,-3,1);
  res_and_msg = lower_bound_check->check();
  REQUIRE(res_and_msg.result==CheckResult::Repairable);

  // repair the field and verify it now passes
  lower_bound_check->repair();
  res_and_msg = lower_bound_check->check();
  REQUIRE(res_and_msg.result==CheckResult::Pass);
}

} // namespace scream
