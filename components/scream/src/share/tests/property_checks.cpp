#include <catch2/catch.hpp>
#include <numeric>

#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_upper_bound_check.hpp"
#include "share/property_checks/field_nan_check.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/grid/point_grid.hpp"

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

  // We don't need a valid grid for these tests
  std::shared_ptr<const AbstractGrid> grid;

  SECTION ("field_not_allocated") {
    Field f(fid);

    // Field must be allocated before being set in the property check
    REQUIRE_THROWS (std::make_shared<FieldNaNCheck>(f,grid));
  }

  SECTION ("field_read_only") {
    Field f(fid);
    f.allocate_view();
    auto cf = f.get_const();

    // Field must not be read-only if repair is needed
    REQUIRE_THROWS (std::make_shared<FieldLowerBoundCheck>(cf,grid,0,true));

    // But ok if no repair is needed
    REQUIRE_NOTHROW (std::make_shared<FieldLowerBoundCheck>(cf,grid,0,false));
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

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_lcols = 2;
  const int nlevs = 12;

  // Create a point grid
  std::shared_ptr<AbstractGrid> grid;
  grid = std::make_shared<PointGrid>("some_grid",num_lcols,nlevs,comm);
  AbstractGrid::dofs_list_type dofs("dogs",grid->get_num_local_dofs());
  AbstractGrid::geo_view_type lat("lat",grid->get_num_local_dofs());
  AbstractGrid::geo_view_type lon("lon",grid->get_num_local_dofs());
  auto lat_h = Kokkos::create_mirror_view(lat);
  auto lon_h = Kokkos::create_mirror_view(lon);
  auto dofs_h = Kokkos::create_mirror_view(dofs);
  for (int i=0; i<grid->get_num_local_dofs(); ++i) {
    dofs_h(i) = num_lcols*comm.rank() + i;
    lat_h(i) = i;
    lon_h(i) = -i;
  }
  Kokkos::deep_copy(dofs,dofs_h);
  Kokkos::deep_copy(lat,lat_h);
  Kokkos::deep_copy(lon,lon_h);
  grid->set_dofs(dofs);
  grid->set_geometry_data("lat",lat);
  grid->set_geometry_data("lon",lon);

  // Create a field
  std::vector<FieldTag> tags = {COL, CMP, LEV};
  std::vector<int> dims = {num_lcols, 3, nlevs};
  FieldIdentifier fid ("field_1",{tags,dims}, m/s,"some_grid");
  Field f(fid);
  f.allocate_view();

  // Check that values are not NaN
  SECTION("field_not_nan_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto nan_check = std::make_shared<FieldNaNCheck>(f,grid);

    // Assign  values to the field and make sure it passes our test for NaNs.
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,neg_pdf);
    f.sync_to_dev();
    auto res_and_msg = nan_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Assign a NaN value to the field, make sure it fails the check,
    auto f_view = f.get_view<Real***,Host>();
    f_view(1,2,3) = std::numeric_limits<Real>::quiet_NaN();
    f.sync_to_dev();
    res_and_msg = nan_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Fail);
    std::string expected_msg =
      "FieldNaNCheck failed.\n"
      "  - field id: " + fid.get_id_string() + "\n"
      "  - entry (1,2,3)\n"
      "  - lat/lon: (1.000000, -1.000000)\n";
    REQUIRE( res_and_msg.msg == expected_msg );
  }

  // Check that the values of a field lie within an interval.
  SECTION ("field_within_interval_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto interval_check = std::make_shared<FieldWithinIntervalCheck>(f, grid, 0, 1, true);
    REQUIRE(interval_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the within-interval check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,pos_pdf);
    f.sync_to_dev();
    auto res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    f.deep_copy(0.5);
    auto f_view = f.get_view<Real***,Host>();
    f_view(1,2,3) = 2.0;
    f_view(0,1,2) = 0.0;
    f.sync_to_dev();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);
    std::string expected_msg =
      "Check failed.\n"
      "  - check name: field_1 within interval [0, 1]\n"
      "  - field id: " + fid.get_id_string() + "\n"
      "  - minimum:\n"
      "    - value: 0\n"
      "    - entry: (0,1,2)\n"
      "    - lat/lon: (0, 0)\n"
      "  - maximum:\n"
      "    - value: 2\n"
      "    - entry: (1,2,3)\n"
      "    - lat/lon: (1, -1)\n";

    REQUIRE(res_and_msg.msg == expected_msg);

    interval_check->repair();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);
  }

  SECTION ("field_within_interval_check_repairable_bounds") {
    using FWIC = FieldWithinIntervalCheck;

    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Repairable bounds must be less tight than bounds
    REQUIRE_THROWS (FWIC(f,grid,0,1,true,1,2));
    REQUIRE_THROWS (FWIC(f,grid,0,1,true,-1,0));

    auto interval_check = std::make_shared<FWIC>(f, grid, 0, 1, true,-1,2);
    REQUIRE(interval_check->can_repair());

    // Assign slightly out-of-bound values to the field
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    ekat::genRandArray(f_data,num_reals,engine,pos_pdf);
    f_data[0] = -0.5;
    f_data[num_reals-1] = 1.5;
    f.sync_to_dev();
    auto res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);

    // Assign completely out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    f_data[0] = -2;
    f_data[num_reals-1] = 3;
    f.sync_to_dev();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Fail);
    std::string expected_msg =
      "Check failed.\n"
      "  - check name: field_1 within interval [0, 1]\n"
      "  - field id: " + fid.get_id_string() + "\n"
      "  - minimum:\n"
      "    - value: -2\n"
      "    - entry: (0,0,0)\n"
      "    - lat/lon: (0, 0)\n"
      "  - maximum:\n"
      "    - value: 3\n"
      "    - entry: (1,2,11)\n"
      "    - lat/lon: (1, -1)\n";

    REQUIRE(res_and_msg.msg == expected_msg);

    interval_check->repair();
    res_and_msg = interval_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);
  }

  // Check that the values of a field are above a lower bound
  SECTION ("field_lower_bound_check") {
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    auto lower_bound_check = std::make_shared<FieldLowerBoundCheck>(f,grid,-1.0, true);
    REQUIRE(lower_bound_check->can_repair());

    // Assign in-bound values to the field and make sure it passes the lower_bound check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = std::numeric_limits<Real>::max() - i*1.0;
    }
    f.sync_to_dev();
    auto res_and_msg = lower_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = -2.0*(i+1);
    }
    f.sync_to_dev();
    res_and_msg = lower_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);
    lower_bound_check->repair();
    res_and_msg = lower_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);
    // Should have repaired to the lower bound:
    f.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f_data[i] == -1.0);
    }
  }

  // Check that the values of a field are above below an upper bound
  SECTION ("field_upper_bound_check") {
    auto upper_bound_check = std::make_shared<FieldUpperBoundCheck>(f,grid,1.0, true);
    REQUIRE(upper_bound_check->can_repair());
    const auto num_reals = f.get_header().get_alloc_properties().get_num_scalars();

    // Assign in-bound values to the field and make sure it passes the upper_bound check
    auto f_data = reinterpret_cast<Real*>(f.get_internal_view_data<Real,Host>());
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = -std::numeric_limits<Real>::max() + i*1.0;
    }
    f.sync_to_dev();
    auto res_and_msg = upper_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);

    // Assign out-of-bounds values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i<num_reals; ++i) {
      f_data[i] = 2.0*(i+1);
    }
    f.sync_to_dev();
    res_and_msg = upper_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Repairable);
    upper_bound_check->repair();
    res_and_msg = upper_bound_check->check();
    REQUIRE(res_and_msg.result==CheckResult::Pass);
    // Should have repaired to the upper bound:
    f.sync_to_host();
    for (int i=0; i<num_reals; ++i) {
      REQUIRE(f_data[i] == 1.0);
    }
  }
}

} // anonymous namespace
