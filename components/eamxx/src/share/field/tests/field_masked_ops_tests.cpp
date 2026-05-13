#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace {

// Helper: build a column mask where even cols -> 1, odd cols -> 0
// using compute_mask on a real field
scream::Field make_col_mask (const scream::Field& f_real, const scream::Field& f_int_tmpl,
                             int ncol)
{
  auto f_tmp = f_real.clone();
  for (int icol = 0; icol < ncol; ++icol) {
    f_tmp.subfield(0, icol).deep_copy(icol % 2 == 0 ? 1 : -1);
  }
  auto mask = f_int_tmpl.clone("mask");
  scream::compute_mask(f_tmp, 0, scream::Comparison::GT, mask);
  return mask;
}

TEST_CASE ("masked_ops") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  int seed = get_random_test_seed();

  const int ncol = 4;  // use 4 cols so even/odd pattern is clear
  const int ncmp = 3;
  const int nlev = 5;

  std::vector<FieldTag> tags = {COL, CMP, LEV};
  std::vector<int>      dims = {ncol, ncmp, nlev};

  FieldIdentifier fid_r ("fr", {tags, dims}, kg, "some_grid", DataType::RealType);
  FieldIdentifier fid_i ("fi", {tags, dims}, kg, "some_grid", DataType::IntType);

  Field f_real(fid_r);
  Field f_int (fid_i);
  f_real.allocate_view();
  f_int.allocate_view();
  randomize_uniform(f_real, seed++);
  randomize_uniform(f_int,  seed++, 0, 100);

  // Convenience: one column's subfield as a reference shape
  auto col0_r = f_real.subfield(0, 0);
  auto col0_i = f_int.subfield(0, 0);

  SECTION ("masked_deep_copy_scalar") {
    // Set dst=1 on masked cols, dst=2 on un-masked cols using negate_mask
    auto dst = f_real.clone();
    dst.deep_copy(0);

    auto mask = make_col_mask(f_real, f_int, ncol);

    dst.deep_copy(1, mask);           // 1 where mask!=0
    dst.deep_copy(2, mask, true);     // 2 where mask==0  (negated)

    auto one = col0_r.clone("one"); one.deep_copy(1);
    auto two = col0_r.clone("two"); two.deep_copy(2);

    for (int icol = 0; icol < ncol; ++icol) {
      auto dsti = dst.subfield(0, icol);
      REQUIRE(views_are_equal(dsti, icol % 2 == 0 ? one : two));
    }
  }

  SECTION ("masked_deep_copy_field") {
    // Copy src field values into dst only at masked columns
    auto src = f_real.clone();
    auto dst = f_real.clone();
    randomize_uniform(src, seed++);
    dst.deep_copy(0);

    auto mask = make_col_mask(f_real, f_int, ncol);

    dst.deep_copy(src, mask);

    auto zero_col = col0_r.clone("zero"); zero_col.deep_copy(0);
    for (int icol = 0; icol < ncol; ++icol) {
      auto dsti = dst.subfield(0, icol);
      auto srci = src.subfield(0, icol);
      if (icol % 2 == 0) {
        REQUIRE(views_are_equal(dsti, srci));
      } else {
        REQUIRE(views_are_equal(dsti, zero_col));
      }
    }
  }

  SECTION ("masked_deep_copy_int_scalar") {
    // Same as masked_deep_copy_scalar but for integer fields
    auto dst = f_int.clone();
    dst.deep_copy(0);

    auto mask = make_col_mask(f_real, f_int, ncol);

    dst.deep_copy(7, mask);
    dst.deep_copy(3, mask, true);

    auto seven = col0_i.clone("seven"); seven.deep_copy(7);
    auto three = col0_i.clone("three"); three.deep_copy(3);

    for (int icol = 0; icol < ncol; ++icol) {
      auto dsti = dst.subfield(0, icol);
      REQUIRE(views_are_equal(dsti, icol % 2 == 0 ? seven : three));
    }
  }

  SECTION ("masked_update") {
    // f2 = alpha*f_real + beta*f2, only on masked cols; unmasked stay at initial value
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_real.clone();
    f2.deep_copy(1);
    f2.update(f_real, 2, 3, mask);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i   = f2.subfield(0, icol);
      auto f_ref = col0_r.clone("ref");
      if (icol % 2 == 0) {
        // ref = 2*f_real + 3
        f_ref.deep_copy(1);
        f_ref.update(f_real.subfield(0, icol), 2, 3);
      } else {
        f_ref.deep_copy(1);
      }
      REQUIRE(views_are_equal(f2i, f_ref));
    }
  }

  SECTION ("masked_update_int") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_int.clone();
    f2.deep_copy(1);
    f2.update(f_int, 2, 3, mask);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i   = f2.subfield(0, icol);
      auto f_ref = col0_i.clone("ref");
      if (icol % 2 == 0) {
        f_ref.deep_copy(1);
        f_ref.update(f_int.subfield(0, icol), 2, 3);
      } else {
        f_ref.deep_copy(1);
      }
      REQUIRE(views_are_equal(f2i, f_ref));
    }
  }

  SECTION ("masked_scale_scalar") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_real.clone();   // starts as a copy of f_real
    f2.scale(2, mask);        // only masked cols multiplied by 2

    auto f_ref = f_real.clone();
    f_ref.scale(2);           // reference: 2*f_real everywhere

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i     = f2.subfield(0, icol);
      auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol)
                                      : f_real.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, f_ref_i));
    }
  }

  SECTION ("masked_scale_field") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2     = f_real.clone();
    auto f_twos = f_real.clone();
    f_twos.deep_copy(2);
    f2.scale(f_twos, mask);

    auto f_ref = f_real.clone();
    f_ref.scale(2);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i     = f2.subfield(0, icol);
      auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol)
                                      : f_real.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, f_ref_i));
    }
  }

  SECTION ("masked_scale_int") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    // Integer scale by field of twos, only on masked cols
    auto f2     = f_int.clone();
    auto f_twos = f_int.clone();
    f_twos.deep_copy(2);
    f2.scale(f_twos, mask);

    auto f_ref = f_int.clone();
    f_ref.scale(f_twos);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i     = f2.subfield(0, icol);
      auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol)
                                      : f_int.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, f_ref_i));
    }
  }

  SECTION ("masked_scale_inv") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2     = f_real.clone();
    auto f_twos = f_real.clone();
    f_twos.deep_copy(2);
    f2.scale_inv(f_twos, mask);

    auto f_ref = f_real.clone();
    f_ref.scale(0.5);   // f_real/2

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i     = f2.subfield(0, icol);
      auto f_ref_i = (icol % 2 == 0) ? f_ref.subfield(0, icol)
                                      : f_real.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, f_ref_i));
    }
  }

  SECTION ("masked_max") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_real.clone(); f2.deep_copy(1);
    auto f3 = f_real.clone(); f3.deep_copy(2);
    f2.max(f3, mask);   // masked cols: max(1,2)=2; unmasked: stay 1

    auto one = col0_r.clone("one"); one.deep_copy(1);
    auto two = col0_r.clone("two"); two.deep_copy(2);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, icol % 2 == 0 ? two : one));
    }
  }

  SECTION ("masked_max_int") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_int.clone(); f2.deep_copy(1);
    auto f3 = f_int.clone(); f3.deep_copy(2);
    f2.max(f3, mask);

    auto one = col0_i.clone("one"); one.deep_copy(1);
    auto two = col0_i.clone("two"); two.deep_copy(2);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, icol % 2 == 0 ? two : one));
    }
  }

  SECTION ("masked_min") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_real.clone(); f2.deep_copy(3);
    auto f3 = f_real.clone(); f3.deep_copy(1);
    f2.min(f3, mask);   // masked cols: min(3,1)=1; unmasked: stay 3

    auto one   = col0_r.clone("one");   one.deep_copy(1);
    auto three = col0_r.clone("three"); three.deep_copy(3);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, icol % 2 == 0 ? one : three));
    }
  }

  SECTION ("masked_min_int") {
    auto mask = make_col_mask(f_real, f_int, ncol);

    auto f2 = f_int.clone(); f2.deep_copy(3);
    auto f3 = f_int.clone(); f3.deep_copy(1);
    f2.min(f3, mask);

    auto one   = col0_i.clone("one");   one.deep_copy(1);
    auto three = col0_i.clone("three"); three.deep_copy(3);

    for (int icol = 0; icol < ncol; ++icol) {
      auto f2i = f2.subfield(0, icol);
      REQUIRE(views_are_equal(f2i, icol % 2 == 0 ? one : three));
    }
  }

  SECTION ("all_masked") {
    // Mask is all 1s: every entry should be updated
    auto mask = f_int.clone("all_ones_mask");
    mask.deep_copy(1);

    auto f2 = f_real.clone(); f2.deep_copy(0);
    auto f3 = f_real.clone(); f3.deep_copy(5);
    f2.deep_copy(5, mask);
    REQUIRE(views_are_equal(f2, f3));

    auto f4 = f_real.clone(); f4.deep_copy(1);
    f4.scale(3, mask);
    auto f_ref = f_real.clone(); f_ref.deep_copy(3);
    REQUIRE(views_are_equal(f4, f_ref));

    auto f5 = f_real.clone(); f5.deep_copy(2);
    auto f6 = f_real.clone(); f6.deep_copy(4);
    f5.max(f6, mask);
    REQUIRE(views_are_equal(f5, f6));
  }

  SECTION ("none_masked") {
    // Mask is all 0s: no entry should be updated
    auto mask = f_int.clone("all_zeros_mask");
    mask.deep_copy(0);

    // deep_copy with zero mask: dst should remain unchanged
    auto dst = f_real.clone(); dst.deep_copy(99);
    dst.deep_copy(1, mask);
    auto ref = f_real.clone(); ref.deep_copy(99);
    REQUIRE(views_are_equal(dst, ref));

    // scale with zero mask: f should remain unchanged
    auto f2 = f_real.clone();
    f2.scale(10, mask);
    REQUIRE(views_are_equal(f2, f_real));

    // max with zero mask: f should remain unchanged
    auto f3 = f_real.clone(); f3.deep_copy(1);
    auto f4 = f_real.clone(); f4.deep_copy(100);
    f3.max(f4, mask);
    auto one = f_real.clone(); one.deep_copy(1);
    REQUIRE(views_are_equal(f3, one));
  }

  SECTION ("nonbinary_mask") {
    // Mask entries !0 should behave as true, regardless of exact value
    auto mask = f_int.clone("nonbinary_mask");
    mask.deep_copy(0);
    // Set even cols to 5 (instead of 1) - should still behave as "masked"
    for (int icol = 0; icol < ncol; ++icol) {
      mask.subfield(0, icol).deep_copy(icol % 2 == 0 ? 5 : 0);
    }

    auto dst = f_real.clone(); dst.deep_copy(0);
    dst.deep_copy(7, mask);

    auto zero = col0_r.clone("zero"); zero.deep_copy(0);
    auto seven = col0_r.clone("seven"); seven.deep_copy(7);

    for (int icol = 0; icol < ncol; ++icol) {
      auto dsti = dst.subfield(0, icol);
      REQUIRE(views_are_equal(dsti, icol % 2 == 0 ? seven : zero));
    }
  }
}

} // anonymous namespace
