#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace {

TEST_CASE ("update") {
  using namespace scream;
  using namespace ekat::units;

  using namespace ShortFieldTagsNames;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  int seed = get_random_test_seed();

  const int ncol = 2;
  const int ncmp = 3;
  const int nlev = 4;

  // Create field
  std::vector<FieldTag> tags = {COL, CMP, LEV};
  std::vector<int>      dims = {ncol,ncmp,nlev};

  FieldIdentifier fid_r ("fr", {tags,dims}, kg, "some_grid", DataType::RealType);
  FieldIdentifier fid_i ("fi", {tags,dims}, kg, "some_grid", DataType::IntType);
  Field f_real (fid_r);
  Field f_int  (fid_i);
  f_real.allocate_view();
  f_int.allocate_view();
  randomize_uniform (f_real,seed++);
  randomize_uniform (f_int, seed++, 0, 100);

  SECTION ("data_type_checks") {
    Field f2 = f_int.clone();

    // Coeffs have wrong data type (precision loss casting to field's data type)
    REQUIRE_THROWS (f2.update (f_int,1.0,1.0));

    // RHS has wrong data type
    REQUIRE_THROWS (f2.update(f_real,1,0));
  }

  SECTION ("deep_copy") {
    SECTION ("real") {
      Field f2 (fid_r);
      f2.allocate_view();

      // Replace f2's content with f_real's content
      f2.deep_copy(f_real);
      REQUIRE (views_are_equal(f2,f_real));
    }
    SECTION ("int") {
      Field f2 (fid_i);
      f2.allocate_view();

      // Replace f2's content with f_int's content
      f2.deep_copy(f_int);
      REQUIRE (views_are_equal(f2,f_int));
    }
  }

  SECTION ("masked_deep_copy") {
    auto f1 = f_real.clone();
    auto f2 = f_real.clone();
    auto f3 = f_real.clone();
    f3.deep_copy(0);
    for (int icol=0; icol<ncol; ++ icol) {
      auto val = icol % 2 == 0 ? 1 : -1;
      f1.subfield(0,icol).deep_copy(val);
    }

    // Compute mask where f1>0 (should be all even cols)
    auto mask = f_int.clone("valid_mask");
    compute_mask(f1,0,Comparison::GT,mask);

    // Set f3=1 where mask=1
    f3.deep_copy(1,mask);
    f3.deep_copy(2,mask,true);

    auto one = f1.subfield(0,0).clone("one");
    auto two = f1.subfield(0,0).clone("two");
    one.deep_copy(1);
    two.deep_copy(2);

    // Check
    for (int icol=0; icol<ncol; ++ icol) {
      auto f3i = f3.subfield(0,icol);
      if (icol % 2 == 0) {
        REQUIRE (views_are_equal(f3i,one));
      } else {
        REQUIRE (views_are_equal(f3i,two));
      }
    }
  }

  SECTION ("scale") {
    SECTION ("real") {
      Field f1 = f_real.clone();
      Field f2 = f_real.clone();

      // x=2, x*y = 2*y
      f1.deep_copy(2.0);
      f1.scale(f2);
      f2.scale(2.0);
      REQUIRE (views_are_equal(f1, f2));
    }

    SECTION ("int") {
      Field f1 = f_int.clone();
      f1.deep_copy(4);
      Field f2 = f_int.clone();
      f2.deep_copy(2);
      Field f3 = f_int.clone();
      f3.deep_copy(2);

      f2.scale(f3);
      REQUIRE (views_are_equal(f1, f2));
    }
  }

  // Helper to build a column mask: even cols -> 1, odd cols -> 0
  auto make_col_mask = [&]() {
    auto f_tmp = f_real.clone();
    for (int icol = 0; icol < ncol; ++icol) {
      f_tmp.subfield(0, icol).deep_copy(icol % 2 == 0 ? 1.0 : -1.0);
    }
    auto mask = f_int.clone("mask");
    compute_mask(f_tmp, 0.0, Comparison::GT, mask);
    return mask;
  };

  SECTION ("masked_update") {
    auto mask = make_col_mask();

    // f2 starts at 1; update(f_real, 2, 3, mask) -> f2 = 2*f_real+3 on masked cols, 1 elsewhere
    auto f2 = f_real.clone();
    f2.deep_copy(1);
    f2.update(f_real, 2, 3, mask);

    auto tgt_col = f_real.subfield(0, 0).clone("tgt");
    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      if (icol % 2 == 0) {
        tgt_col.deep_copy(f2i);
        tgt_col.update(f_real.subfield(0,icol),2,3);
      } else {
        tgt_col.deep_copy(0);
      }
      REQUIRE (views_are_equal(f2i, tgt_col));
    }
  }

  SECTION ("masked_scale") {
    auto mask = make_col_mask();

    SECTION ("scalar") {
      // scale by 2 on masked cols only; unmasked cols unchanged
      auto f2 = f_real.clone();
      f2.scale(2.0, mask);

      auto f_ref = f_real.clone();
      f_ref.scale(2.0);  // 2*f_real everywhere (reference for masked cols)
      for (int icol = 0; icol < ncol; ++icol) {
        auto f2i   = f_real.subfield(0, icol).clone();
        auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol) : f_real.subfield(0, icol);
        REQUIRE (views_are_equal(f2i, f_ref_i));
      }
    }

    SECTION ("field") {
      // scale by a field of twos on masked cols only
      auto f2 = f_real.clone();
      auto f3 = f_real.clone();
      f3.deep_copy(2.0);
      f2.scale(f3, mask);

      auto f_ref = f_real.clone();
      f_ref.scale(2.0);
      for (int icol = 0; icol < ncol; ++icol) {
        auto f2i   = f2.subfield(0, icol);
        auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol) : f_real.subfield(0, icol);
        REQUIRE (views_are_equal(f2i, f_ref_i));
      }
    }
  }

  SECTION ("masked_scale_inv") {
    auto mask = make_col_mask();

    // scale_inv by 2 on masked cols: f2 = f_real/2 on masked, f_real elsewhere
    auto f2 = f_real.clone();
    auto f3 = f_real.clone();
    f3.deep_copy(2.0);
    f2.scale_inv(f3, mask);

    auto f_ref = f_real.clone();
    f_ref.scale(0.5);
    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i   = f2.subfield(0, icol);
      auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol) : f_real.subfield(0, icol);
      REQUIRE (views_are_equal(f2i, f_ref_i));
    }
  }

  SECTION ("masked_max") {
    auto mask = make_col_mask();

    // f2=1, f3=2; max(f2,f3,mask) -> masked cols become 2, unmasked stay 1
    auto f2 = f_real.clone();
    auto f3 = f_real.clone();
    f2.deep_copy(1.0);
    f3.deep_copy(2.0);
    f2.max(f3, mask);

    auto one = f_real.subfield(0, 0).clone("one");
    auto two = f_real.subfield(0, 0).clone("two");
    one.deep_copy(1.0);
    two.deep_copy(2.0);
    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE (views_are_equal(f2i, icol % 2 == 0 ? two : one));
    }
  }

  SECTION ("masked_min") {
    auto mask = make_col_mask();

    // f2=3, f3=1; min(f2,f3,mask) -> masked cols become 1, unmasked stay 3
    auto f2 = f_real.clone();
    auto f3 = f_real.clone();
    f2.deep_copy(3.0);
    f3.deep_copy(1.0);
    f2.min(f3, mask);

    auto one   = f_real.subfield(0, 0).clone("one");
    auto three = f_real.subfield(0, 0).clone("three");
    one.deep_copy(1.0);
    three.deep_copy(3.0);
    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE (views_are_equal(f2i, icol % 2 == 0 ? one : three));
    }
  }

  SECTION ("max-min") {
    SECTION ("real") {
      Field one = f_real.clone();
      Field two = f_real.clone();
      one.deep_copy(1.0);
      two.deep_copy(2.0);

      Field f1 = one.clone();
      Field f2 = two.clone();
      f1.max(f2);
      REQUIRE (views_are_equal(f1, f2));

      Field f3 = one.clone();
      Field f4 = two.clone();
      f4.min(f3);
      REQUIRE (views_are_equal(f3, f4));

      // Check that updating with rhs==fill_value ignores the rhs
      f3.deep_copy(constants::fill_value<Real>);
      f3.get_header().set_may_be_filled(true);
      f2.deep_copy(1.0);
      f2.max(f3);
      REQUIRE (views_are_equal(f2,one));
    }

    SECTION ("int") {
      Field one = f_int.clone();
      Field two = f_int.clone();
      one.deep_copy(1);
      two.deep_copy(2);

      Field f1 = one.clone();
      Field f2 = two.clone();
      f1.max(f2);
      REQUIRE (views_are_equal(f1, f2));

      Field f3 = one.clone();
      Field f4 = two.clone();
      f4.min(f3);
      REQUIRE (views_are_equal(f3, f4));

      // Check that updating with rhs==fill_value ignores the rhs
      f3.deep_copy(constants::fill_value<int>);
      f3.get_header().set_may_be_filled(true);
      f2.deep_copy(1);
      f2.max(f3);
      REQUIRE (views_are_equal(f2,one));
    }
  }

  SECTION ("scale_inv") {
    SECTION ("real") {
      Field f1 = f_real.clone();
      Field f2 = f_real.clone();
      Field f3 = f_real.clone();

      f3.deep_copy(2.0);
      f1.scale(f3);
      f3.deep_copy(0.5);
      f2.scale_inv(f3);
      REQUIRE (views_are_equal(f1, f2));
    }

    SECTION ("int") {
      Field f1 = f_int.clone();
      f1.deep_copy(4);
      Field f2 = f_int.clone();
      f2.deep_copy(2);

      f1.scale_inv(f2);
      REQUIRE (views_are_equal(f1, f2));
    }
  }

  SECTION ("update") {
    SECTION ("real") {
      Field f2 = f_real.clone();
      Field f3 = f_real.clone();

      // x+x == 2*x
      f2.update(f_real,1,1);
      f3.scale(2);
      REQUIRE (views_are_equal(f2,f3));

      // Adding 2*f_real to N*f3 should give 2*f_real (f3==0)
      f3.deep_copy(0.0);
      f3.update(f_real,2,10);
      REQUIRE (views_are_equal(f3,f2));

      // Same, but we discard current content of f3
      f3.update(f_real,2,0);
      REQUIRE (views_are_equal(f3,f2));

      // Check that updating with rhs==fill_value ignores the rhs
      Field one = f_real.clone();
      one.deep_copy(1.0);

      f3.deep_copy(constants::fill_value<Real>);
      f3.get_header().set_may_be_filled(true);
      f2.deep_copy(1.0);
      f2.update(f3,1,1);
      if (not views_are_equal(f2,one)) {
        print_field_hyperslab(f2);
      }
      REQUIRE (views_are_equal(f2,one));
    }

    SECTION ("int") {
      Field f2 = f_int.clone();
      Field f3 = f_int.clone();

      // x+x == 2*x
      f2.update(f_int,1,1);
      f3.scale(2);
      REQUIRE (views_are_equal(f2,f3));

      // Adding 2*f_int to N*f3 should give 2*f_int (f3==0)
      f3.deep_copy(0);
      f3.update(f_int,2,10);
      REQUIRE (views_are_equal(f3,f2));

      // Same, but we discard current content of f3
      f3.update(f_int,2,0);
      REQUIRE (views_are_equal(f3,f2));

      // Check that updating with rhs==fill_value ignores the rhs
      Field one = f_int.clone();
      one.deep_copy(1);

      f3.deep_copy(constants::fill_value<int>);
      f3.get_header().set_may_be_filled(true);
      f2.deep_copy(1);
      f2.update(f3,1,1);
      REQUIRE (views_are_equal(f2,one));
    }
  }
}

} // anonymous namespace
