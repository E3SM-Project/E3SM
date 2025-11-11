#include <catch2/catch.hpp>

#include "pc_tests_helpers.hpp"

#include "share/property_checks/field_nan_check.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("exceptions") {
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grid
  const int num_lcols = 10;
  const int nlevs = 16;
  auto grid = create_test_grid(comm,num_lcols,nlevs);

  FieldIdentifier fid("f",grid->get_2d_scalar_layout(),ekat::units::m,"");
  SECTION ("field_not_allocated") {
    Field f(fid);

    REQUIRE_THROWS (std::make_shared<FieldNaNCheck>(f,grid));
  }

  SECTION ("field_read_only") {
    Field f(fid);
    f.allocate_view();
    auto cf = f.get_const();

    // Field must not be read-only if repair is needed
    REQUIRE_THROWS (std::make_shared<FieldNaNCheck>(cf,grid));
  }
}

TEST_CASE("nan_check") {
  using namespace ShortFieldTagsNames;

  auto seed = get_random_test_seed();

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_lcols = 2;
  const int nlevs = 12;

  // Create grid
  auto grid = create_test_grid(comm,num_lcols,nlevs);

  // Create field to test, and extra data field
  auto f = create_test_field(grid);
  randomize_uniform(f,seed,0.01,0.99);
  auto data = create_data_field(grid);

  // Create the check
  auto nan_check = std::make_shared<FieldNaNCheck>(f,grid);
  nan_check->set_additional_data_field(data);

  // Run the check with the field being ok
  auto res_and_msg = nan_check->check();
  REQUIRE(res_and_msg.result==CheckResult::Pass);

  // Run the check with the field being NOT ok
  auto nan = std::numeric_limits<Real>::quiet_NaN();
  auto f_123 = f.subfield(COL,1).subfield(CMP,2).subfield(LEV,3);
  f_123.deep_copy(nan);

  res_and_msg = nan_check->check();
  REQUIRE(res_and_msg.result==CheckResult::Fail);

  std::string expected_msg =
    "FieldNaNCheck failed.\n"
    "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "  - indices (w/ global column index): (1,2,3)\n"
    "  - lat/lon: (1.000000, -1.000000)\n"
    "  - additional data (w/ local column index):\n\n"
    "     data<ncol>(2)\n\n"+
    "  data(1)\n"+
    "    1, \n\n"+
    "  END OF ADDITIONAL DATA\n";

  REQUIRE( res_and_msg.msg == expected_msg );
}

} // namespace scream
