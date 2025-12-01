#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

#include <limits>

namespace scream {

TEST_CASE("field_contractions") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // The following two functions are used in both horiz_contraction and
  // vert_contraction below
  auto sum_n    = [](int n) { return n * (n + 1) / 2; };
  auto sum_n_sq = [](int n) { return n * (n + 1) * (2 * n + 1) / 6; };

  auto seed = get_random_test_seed();

  SECTION("horiz_contraction") {
    // A numerical tolerance
    // Accumulations in the Kokkos threaded reductions may be done in a
    // different order than the manual ones below, so we can only test
    // correctness up to a tolerance
    auto tol = std::numeric_limits<Real>::epsilon() * 100;

    int dim0 = 129;
    int dim1 = 4;
    int dim2 = 17;

    // Set a weight field
    FieldIdentifier f00("f", {{COL}, {dim0}}, m / s, "g");
    Field field00(f00);
    field00.allocate_view();
    field00.sync_to_host();
    auto v00 = field00.get_strided_view<Real *, Host>();
    for(int i = 0; i < dim0; ++i) {
      // By design, denominator is the sum of the first dim0 integers
      v00(i) = sp(i + 1) / sp(sum_n(dim0));
    }
    field00.sync_to_dev();

    // Create (random) sample fields
    FieldIdentifier fsc("f", {{}, {}}, m / s, "g");  // scalar
    FieldIdentifier f10("f", {{COL, CMP}, {dim0, dim1}}, m / s, "g");
    FieldIdentifier f11("f", {{COL, LEV}, {dim0, dim2}}, m / s, "g");
    FieldIdentifier f20("f", {{COL, CMP, LEV}, {dim0, dim1, dim2}}, m / s, "g");
    Field fieldsc(fsc);
    Field field10(f10);
    Field field11(f11);
    Field field20(f20);
    fieldsc.allocate_view();
    field10.allocate_view();
    field11.allocate_view();
    field20.allocate_view();
    randomize_uniform(fieldsc, seed++);
    randomize_uniform(field10, seed++);
    randomize_uniform(field11, seed++);
    randomize_uniform(field20, seed++);

    FieldIdentifier F_x("fx", {{COL}, {dim0}}, m / s, "g");
    FieldIdentifier F_y("fy", {{LEV}, {dim2}}, m / s, "g");
    FieldIdentifier F_z("fz", {{CMP}, {dim1}}, m / s, "g");
    FieldIdentifier F_w("fyz", {{CMP, LEV}, {dim1, dim2}}, m / s, "g");

    Field field_x(F_x);
    Field field_y(F_y);
    Field field_z(F_z);
    Field field_w(F_w);

    // Test invalid inputs
    REQUIRE_THROWS(horiz_contraction(fieldsc, field_x, field00));  // x not allocated yet

    field_x.allocate_view();
    field_y.allocate_view();
    field_z.allocate_view();
    field_w.allocate_view();

    REQUIRE_THROWS(horiz_contraction(fieldsc, field_y, field_x));  // unmatching layout
    REQUIRE_THROWS(horiz_contraction(field_z, field11, field11));  // wrong weight layout

    Field result;

    // Ensure a scalar case works
    result = fieldsc.clone();
    horiz_contraction(result, field00, field00);
    result.sync_to_host();
    auto v = result.get_view<Real, Host>();
    // The numerator is the sum of the squares of the first dim0 integers
    // The denominator is the sum of the first dim0 integers squared
    Real wavg = sp(sum_n_sq(dim0)) / sp(sum_n(dim0) * sum_n(dim0));
    REQUIRE_THAT(v(), Catch::Matchers::WithinRel(wavg, tol));

    // Repeat but with masked values
    result = fieldsc.clone();
    // inject a mask as the last entry
    auto field00_masked = field00.clone();
    auto mask_of_field00 = field00_masked.clone();
    mask_of_field00.deep_copy(sp(1.0));
    mask_of_field00.sync_to_host();
    auto mask = mask_of_field00.get_view<Real *, Host>();
    mask(dim0 - 1) = sp(0.0);
    mask_of_field00.sync_to_dev();
    field00_masked.get_header().set_extra_data("mask_data", mask_of_field00);
    field00_masked.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
    field00_masked.sync_to_dev();
    auto result_mask = result.clone();
    result.get_header().set_extra_data("mask_data", result_mask);
    result.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
    horiz_contraction(result, field00_masked, field00);
    result.sync_to_host();
    v = result.get_view<Real, Host>();
    Real wavg_sum1 = 0;
    Real wavg_sum2 = 0;
    auto wavg_v00 = field00.get_view<const Real *, Host>();
    for(int i = 0; i < dim0; ++i) {
      wavg_sum1 += mask(i) * wavg_v00(i) * wavg_v00(i);
      wavg_sum2 += mask(i) * wavg_v00(i);
    }
    REQUIRE_THAT(v(), Catch::Matchers::WithinRel(wavg_sum1/wavg_sum2, tol));

    // Repeat but with ALL masked values
    result = fieldsc.clone();
    // inject a mask as the last entry
    field00_masked = field00.clone();
    mask_of_field00 = field00_masked.clone();
    mask_of_field00.deep_copy(sp(1.0));
    mask_of_field00.sync_to_host();
    mask = mask_of_field00.get_view<Real *, Host>();
    for(int i = 0; i < dim0; ++i) {
      mask(i) = sp(0.0);
    }
    mask_of_field00.sync_to_dev();
    field00_masked.get_header().set_extra_data("mask_data", mask_of_field00);
    Real mask_v = constants::fill_value<Real>;
    field00_masked.get_header().set_extra_data("mask_value", mask_v);
    field00_masked.sync_to_dev();
    result_mask = result.clone();
    result.get_header().set_extra_data("mask_data", result_mask);
    result.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
    horiz_contraction(result, field00_masked, field00);
    result.sync_to_host();
    v = result.get_view<Real, Host>();
    REQUIRE(v() == mask_v);

    // Test higher-order cases
    result = field_z.clone();
    horiz_contraction(result, field10, field00);
    REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
            std::vector<FieldTag>({CMP}));
    REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim1);

    result = field_y.clone();
    horiz_contraction(result, field11, field00);
    REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
            std::vector<FieldTag>({LEV}));
    REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim2);

    result = field_w.clone();
    horiz_contraction(result, field20, field00);
    REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
            std::vector<FieldTag>({CMP, LEV}));
    REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim1);
    REQUIRE(result.get_header().get_identifier().get_layout().dim(1) == dim2);

    // Check a 3D case
    field20.sync_to_host();
    result.sync_to_host();
    auto manual_result = result.clone();
    manual_result.deep_copy(0);
    manual_result.sync_to_host();
    auto v2 = field20.get_strided_view<Real ***, Host>();
    auto mr = manual_result.get_strided_view<Real **, Host>();
    auto rr = result.get_strided_view<Real **, Host>();

    for(int j = 0; j < dim1; ++j) {
      for(int k = 0; k < dim2; ++k) {
        for(int i = 0; i < dim0; ++i) {
          mr(j, k) += v00(i) * v2(i, j, k);
        }
        REQUIRE_THAT(rr(j, k), Catch::Matchers::WithinRel(mr(j, k), tol));
      }
    }
  }

  SECTION("vert_contraction") {
    // A numerical tolerance
    // Accumulations in the Kokkos threaded reductions may be done in a
    // different order than the manual ones below, so we can only test
    // correctness up to a tolerance
    auto tol = std::numeric_limits<Real>::epsilon() * 100;

    std::vector<FieldTag> lev_tags = {LEV, ILEV};
    // iterate over lev_tags
    for(auto lev_tag : lev_tags) {
      int dim0 = 18;
      int dim1 = 9;
      // Note that parallel reduction is happening over dim2 (LEV/ILEV)
      int dim2 = lev_tag == LEV ? 225 : 226;

      // Set a weight field
      FieldIdentifier f00("f", {{lev_tag}, {dim2}}, m / s, "g");
      Field field00(f00);
      field00.allocate_view();
      field00.sync_to_host();
      auto v00 = field00.get_strided_view<Real *, Host>();
      for(int i = 0; i < dim2; ++i) {
        // The denominator is the sum of the first dim2 integers (analytically
        // known)
        v00(i) = sp(i + 1) / sp(sum_n(dim2));
      }
      field00.sync_to_dev();

      // Create (random) sample fields
      FieldIdentifier fsc("f", {{}, {}}, m / s, "g");  // scalar
      FieldIdentifier f10("f", {{COL, lev_tag}, {dim0, dim2}}, m / s, "g");
      FieldIdentifier f11("f", {{CMP, lev_tag}, {dim1, dim2}}, m / s, "g");
      FieldIdentifier f20("f", {{COL, CMP, lev_tag}, {dim0, dim1, dim2}}, m / s,
                          "g");
      Field fieldsc(fsc);
      Field field10(f10);
      Field field11(f11);
      Field field20(f20);
      fieldsc.allocate_view();
      field10.allocate_view();
      field11.allocate_view();
      field20.allocate_view();
      randomize_uniform(fieldsc, seed++);
      randomize_uniform(field10, seed++);
      randomize_uniform(field11, seed++);
      randomize_uniform(field20, seed++);

      FieldIdentifier F_x("fx", {{COL}, {dim0}}, m / s, "g");
      FieldIdentifier F_y("fy", {{CMP}, {dim1}}, m / s, "g");
      FieldIdentifier F_z("fz", {{COL, CMP}, {dim0, dim1}}, m / s, "g");

      Field field_x(F_x);
      Field field_y(F_y);
      Field field_z(F_z);

      // Test invalid inputs
      REQUIRE_THROWS(vert_contraction(fieldsc, field_x, field00));  // x not allocated yet

      field_x.allocate_view();
      field_y.allocate_view();
      field_z.allocate_view();

      REQUIRE_THROWS(vert_contraction(fieldsc, field_y, field_x));  // unmatching layout
      REQUIRE_THROWS(vert_contraction(field_z, field11, field11));  // wrong weight layout

      Field result;

      // Add test for invalid rank-2 weight field layout
      FieldIdentifier bad_w("bad_w", {{CMP, lev_tag}, {dim1, dim2}}, m / s,
                            "g");
      Field bad_weight(bad_w);
      bad_weight.allocate_view();
      REQUIRE_THROWS(vert_contraction(result, field20, bad_weight));

      // Add test for mismatched weight field dimensions
      FieldIdentifier wrong_size_w(
          "wrong_w", {{COL, lev_tag}, {dim0 + 1, dim2}}, m / s, "g");
      Field wrong_weight(wrong_size_w);
      wrong_weight.allocate_view();
      REQUIRE_THROWS(vert_contraction(result, field20, wrong_weight));

      // Ensure a scalar case works
      result = fieldsc.clone();
      vert_contraction(result, field00, field00);
      result.sync_to_host();
      auto v = result.get_view<Real, Host>();
      // The numerator is the sum of the squares of the first dim2 integers
      // (analytically known). The denominator is the sum of the first dim2
      // integers squared (analytically known)
      Real havg = sp(sum_n_sq(dim2)) / sp(sum_n(dim2) * sum_n(dim2));
      REQUIRE_THAT(v(), Catch::Matchers::WithinRel(havg, tol));

      // Repeat but with masked values
      result = fieldsc.clone();
      // inject a mask as the last entry
      auto field00_masked = field00.clone();
      auto mask_of_field00 = field00_masked.clone();
      mask_of_field00.deep_copy(sp(1.0));
      mask_of_field00.sync_to_host();
      auto mask = mask_of_field00.get_view<Real *, Host>();
      mask(dim2 - 1) = sp(0.0);
      mask_of_field00.sync_to_dev();
      field00_masked.get_header().set_extra_data("mask_data", mask_of_field00);
      field00_masked.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
      field00_masked.sync_to_dev();
      auto result_mask = result.clone();
      result.get_header().set_extra_data("mask_data", result_mask);
      result.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
      vert_contraction(result, field00_masked, field00, true);
      result.sync_to_host();
      v = result.get_view<Real, Host>();
      Real wavg_sum1 = sp(0.0);
      Real wavg_sum2 = sp(0.0);
      auto wavg_v00 = field00.get_view<const Real *, Host>();
      for(int i = 0; i < dim2; ++i) {
        wavg_sum1 += mask(i) * wavg_v00(i) * wavg_v00(i);
        wavg_sum2 += mask(i) * wavg_v00(i);
      }
      REQUIRE_THAT(v(), Catch::Matchers::WithinRel(wavg_sum1/wavg_sum2, tol));

      // Repeat but with ALL masked values
      result = fieldsc.clone();
      // inject a mask as the last entry
      field00_masked = field00.clone();
      mask_of_field00 = field00_masked.clone();
      mask_of_field00.deep_copy(sp(1.0));
      mask_of_field00.sync_to_host();
      mask = mask_of_field00.get_view<Real *, Host>();
      for (int i=0; i < dim2; ++i) {
        mask(i) = sp(0.0);
      }
      mask_of_field00.sync_to_dev();
      field00_masked.get_header().set_extra_data("mask_data", mask_of_field00);
      Real mask_v = constants::fill_value<Real>;
      field00_masked.get_header().set_extra_data("mask_value", mask_v);
      field00_masked.sync_to_dev();
      result_mask = result.clone();
      result.get_header().set_extra_data("mask_data", result_mask);
      result.get_header().set_extra_data("mask_value", constants::fill_value<Real>);
      vert_contraction(result, field00_masked, field00, true);
      result.sync_to_host();
      v = result.get_view<Real, Host>();
      REQUIRE(v() == mask_v);

      // Test higher-order cases
      result = field_x.clone();
      vert_contraction(result, field10, field00);
      REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
              std::vector<FieldTag>({COL}));
      REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim0);

      // Check a 2D case with 1D weight
      field10.sync_to_host();
      result.sync_to_host();
      auto manual_result = result.clone();
      manual_result.deep_copy(0);
      manual_result.sync_to_host();
      auto v1 = field10.get_strided_view<Real **, Host>();
      auto mr = manual_result.get_strided_view<Real *, Host>();
      auto rr = result.get_strided_view<Real *, Host>();
      for(int i = 0; i < dim0; ++i) {
        for(int j = 0; j < dim2; ++j) {
          mr(i) += v00(j) * v1(i, j);
        }
        REQUIRE_THAT(rr(i), Catch::Matchers::WithinRel(mr(i), tol));
      }

      result = field_y.clone();
      vert_contraction(result, field11, field00);
      REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
              std::vector<FieldTag>({CMP}));
      REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim1);

      result = field_z.clone();
      vert_contraction(result, field20, field00);
      REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
              std::vector<FieldTag>({COL, CMP}));
      REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim0);
      REQUIRE(result.get_header().get_identifier().get_layout().dim(1) == dim1);

      // Check a 3D case with 1D weight
      field20.sync_to_host();
      result.sync_to_host();
      manual_result = result.clone();
      manual_result.deep_copy(0);
      manual_result.sync_to_host();
      auto v2  = field20.get_strided_view<Real ***, Host>();
      auto mr2 = manual_result.get_strided_view<Real **, Host>();
      auto rr2 = result.get_strided_view<Real **, Host>();
      for(int i = 0; i < dim0; ++i) {
        for(int j = 0; j < dim1; ++j) {
          for(int k = 0; k < dim2; ++k) {
            mr2(i, j) += v00(k) * v2(i, j, k);
          }
          REQUIRE_THAT(rr2(i, j), Catch::Matchers::WithinRel(mr2(i, j), tol));
        }
      }

      // Check a 3D case with 2D weight
      result = field_z.clone();
      vert_contraction(result, field20, field10);
      REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
              std::vector<FieldTag>({COL, CMP}));
      REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim0);
      REQUIRE(result.get_header().get_identifier().get_layout().dim(1) == dim1);

      field20.sync_to_host();
      result.sync_to_host();
      manual_result = result.clone();
      manual_result.deep_copy(0);
      manual_result.sync_to_host();
      auto mr3 = manual_result.get_strided_view<Real **, Host>();
      auto rr3 = result.get_strided_view<Real **, Host>();
      for(int i = 0; i < dim0; ++i) {
        for(int j = 0; j < dim1; ++j) {
          for(int k = 0; k < dim2; ++k) {
            mr3(i, j) += v1(i, k) * v2(i, j, k);
          }
          REQUIRE_THAT(rr3(i, j), Catch::Matchers::WithinRel(mr3(i, j), tol));
        }
      }
    }
  }
}

} // namespace scream
