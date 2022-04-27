#include <catch2/catch.hpp>
#include <numeric>

#include "share/property_checks/field_positivity_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_upper_bound_check.hpp"
#include "share/property_checks/field_nan_check.hpp"
#include "share/property_checks/check_and_repair_wrapper.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("property_check_base", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags = {EL, GP, LEV};
  std::vector<int> dims = {2, 3, 12};
  FieldIdentifier fid ("field_1",{tags,dims}, m/s,"some_grid");

  SECTION ("field_not_allocated") {
    Field f(fid);

    // Field must be allocated before being set in the property check
    REQUIRE_THROWS (std::make_shared<FieldNaNCheck>(f));
  }

  SECTION ("field_read_only") {
    Field f(fid);
    f.allocate_view();
    auto cf = f.get_const();

    // Field must not be read-only if repair is needed
    REQUIRE_THROWS (std::make_shared<FieldPositivityCheck>(cf,true));

    // But ok if no repair is needed
    REQUIRE_NOTHROW (std::make_shared<FieldPositivityCheck>(cf,false));
  }
}

TEST_CASE("property_checks", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto engine = setup_random_test();
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pos_pdf(0.01,0.99);
  RPDF neg_pdf(-0.99, -0.01);

  std::vector<FieldTag> tags = {EL, GP, LEV};
  std::vector<int> dims = {2, 3, 12};
  FieldIdentifier fid ("field_1",{tags,dims}, m/s,"some_grid");
  Field f(fid);
  f.allocate_view();

  // Check positivity.
  SECTION ("field_positivity_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Assign positive values to the field and make sure it passes our test for
    // positivity.
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,pos_pdf);
    f.sync_to_dev();

    auto positivity_check = std::make_shared<FieldPositivityCheck>(f,false);
    REQUIRE(not positivity_check->can_repair());
    auto checkResult = positivity_check->check();
    REQUIRE(checkResult.pass);

    // Assign non-positive values to the field and make sure it fails the check.
    ekat::genRandArray(f_data,num_reals,engine,neg_pdf);
    f.sync_to_dev();
    checkResult = positivity_check->check();
    REQUIRE(not checkResult.pass);
  }

  // Check positivity with repairs.
  SECTION ("field_positivity_check_with_repairs") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Assign non-positive values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,neg_pdf);
    f.sync_to_dev();

    auto positivity_check = std::make_shared<FieldPositivityCheck>(f,true);
    REQUIRE(positivity_check->can_repair());
    auto checkResult = positivity_check->check();
    REQUIRE(not checkResult.pass);
    positivity_check->repair();
    checkResult = positivity_check->check();
    REQUIRE(checkResult.pass);
  }

  // Check that values are not NaN
  SECTION("field_not_nan_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto nan_check = std::make_shared<FieldNaNCheck>(f);

    // Assign  values to the field and make sure it passes our test for NaNs.
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,neg_pdf);
    f.sync_to_dev();
    auto checkResult = nan_check->check();
    REQUIRE(checkResult.pass);

    // Assign a NaN value to the field, make sure it fails the check,
    auto f_view = f.get_view<Real***,Host>();
    f_view(1,2,3) = std::numeric_limits<Real>::quiet_NaN();
    f.sync_to_dev();
    checkResult = nan_check->check();
    REQUIRE(not checkResult.pass);
    std::string expected_msg = 
      "FieldNaNCheck failed.\n"
      "  Invalid values found at position (1,2,3).\n";

    REQUIRE( checkResult.msg == expected_msg );
  }

  // Check that the values of a field lie within an interval.
  SECTION ("field_within_interval_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto interval_check = std::make_shared<FieldWithinIntervalCheck>(f, 0, 1, true);
    REQUIRE(interval_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the within-interval check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,pos_pdf);
    f.sync_to_dev();
    auto checkResult = interval_check->check();
    REQUIRE(checkResult.pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] *= -1;
    }
    f.sync_to_dev();
    checkResult = interval_check->check();
    REQUIRE(not checkResult.pass);
    interval_check->repair();
    checkResult = interval_check->check();
    REQUIRE(checkResult.pass);
  }

  // Check that the values of a field are above a lower bound
  SECTION ("field_lower_bound_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto lower_bound_check = std::make_shared<FieldLowerBoundCheck>(f,-1.0, true);
    REQUIRE(lower_bound_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the lower_bound check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = std::numeric_limits<Real>::max() - i*1.0; 
    }
    f.sync_to_dev();
    auto checkResult = lower_bound_check->check();
    REQUIRE(checkResult.pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = -2.0*(i+1);
    }
    f.sync_to_dev();
    checkResult = lower_bound_check->check();
    REQUIRE(not checkResult.pass);
    lower_bound_check->repair();
    checkResult = lower_bound_check->check();
    REQUIRE(checkResult.pass);
    // Should have repaired to the lower bound:
    f.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f_data[i] == -1.0);
    }
  }

  // Check that the values of a field are above below an upper bound
  SECTION ("field_upper_bound_check") {
    auto upper_bound_check = std::make_shared<FieldUpperBoundCheck>(f,1.0, true);
    REQUIRE(upper_bound_check->can_repair());
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Assign in-bound values to the field and make sure it passes the upper_bound check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = -std::numeric_limits<Real>::max() + i*1.0; 
    }
    f.sync_to_dev();
    auto checkResult = upper_bound_check->check();
    REQUIRE(checkResult.pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = 2.0*(i+1);
    }
    f.sync_to_dev();
    checkResult = upper_bound_check->check();
    REQUIRE(not checkResult.pass);
    upper_bound_check->repair();
    checkResult = upper_bound_check->check();
    REQUIRE(checkResult.pass);
    // Should have repaired to the upper bound:
    f.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f_data[i] == 1.0);
    }
  }

  SECTION ("check_and_repair_wrapper") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Two separate FPC for check and for repair
    constexpr Real ub_check  = 1.0;
    constexpr Real ub_repair = 0.0;
    auto check  = std::make_shared<FieldUpperBoundCheck>(f,ub_check);
    auto repair = std::make_shared<FieldUpperBoundCheck>(f,ub_repair,true);

    auto check_and_repair = std::make_shared<CheckAndRepairWrapper>(check, repair);
    REQUIRE(check_and_repair->can_repair());

    // Assign out-of-bound values to the field, and ensure check fails
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = 2.0;
    }
    f.sync_to_dev();
    auto checkResult = check_and_repair->check();
    REQUIRE(not checkResult.pass);
    // Repair the field, and make sure the field values match ub_repair
    check_and_repair->repair();
    f.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f_data[i] == ub_repair);
    }
  }
}

} // anonymous namespace
