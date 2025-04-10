#include <catch2/catch.hpp>
#include <numeric>

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "share/grid/point_grid.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("utils") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;
  using kt_host = KokkosTypes<HostDevice>;

  using P8 = ekat::Pack<Real,8>;

  ekat::Comm comm(MPI_COMM_WORLD);
  bool parallel_test = comm.size()>1;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");
  Field f1(fid);
  f1.get_header().get_alloc_properties().request_allocation(P8::n);
  f1.allocate_view();

  SECTION("compare-rank-0") {
    // create two fields with rank-0 and get their views
    std::vector<FieldTag> tags_0 = {};
    std::vector<int> dims_0      = {};
    FieldIdentifier fid_01("field_01", {tags_0, dims_0}, m / s, "some_grid");
    FieldIdentifier fid_02("field_02", {tags_0, dims_0}, m / s, "some_grid");
    Field f01(fid_01);
    Field f02(fid_02);
    f01.allocate_view();
    f02.allocate_view();
    auto f01v = f01.get_view<Real>();
    auto f02v = f02.get_view<Real>();
    // fill the views with the same values
    Real val = 54321;
    Kokkos::deep_copy(f01v, val);
    Kokkos::deep_copy(f02v, val);
    // check that the views are equal
    REQUIRE(views_are_equal(f01, f02));

    // fill the views with different values
    Kokkos::deep_copy(f02v, 1 / val);
    // check that the views are not equal
    REQUIRE(not views_are_equal(f01, f02));
  }

  SECTION ("compare") {

    Field f2(fid);
    f2.allocate_view();

    auto v1 = f1.get_strided_view<P8**>();
    auto v2 = f2.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto am_i_root = comm.am_i_root();
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v2(i,j) = i*dim1+j;

      int jpack = j / P8::n;
      int jvec = j % P8::n;
      v1(i,jpack)[jvec] = i*dim1+j;
      if (parallel_test and not am_i_root) {
        v1(i,jpack)[jvec] *= -1;
      }
    });
    Kokkos::fence();

    if (not parallel_test) {
      // The views were filled the same way, so they should test equal
      REQUIRE(views_are_equal(f1,f2));
    } else {
      // On rank 0, all is equal, but on other ranks all is different
      if (comm.am_i_root()) {
        REQUIRE(views_are_equal(f1,f2));
      } else {
        REQUIRE(not views_are_equal(f1,f2));
      }

      // On all ranks, the *global* comparison should fail
      REQUIRE(not views_are_equal(f1,f2,&comm));
    }
  }

  SECTION ("sum") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    Real lsum = offset*lsize + lsize*(lsize+1) / 2.0;
    Real gsum = gsize*(gsize+1) / 2.0;

    REQUIRE(field_sum<Real>(f1)==lsum);
    REQUIRE(field_sum<Real>(f1,&comm)==gsum);
  }

  // The following two functions are used in both horiz_contraction and
  // vert_contraction below
  auto sum_n    = [](int n) { return n * (n + 1) / 2; };
  auto sum_n_sq = [](int n) { return n * (n + 1) * (2 * n + 1) / 6; };

  SECTION("horiz_contraction") {
    // A numerical tolerance
    // Accumulations in the Kokkos threaded reductions may be done in a
    // different order than the manual ones below, so we can only test
    // correctness up to a tolerance
    auto tol = std::numeric_limits<Real>::epsilon() * 100;

    using RPDF  = std::uniform_real_distribution<Real>;
    auto engine = setup_random_test();
    RPDF pdf(0, 1);

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
    randomize(fieldsc, engine, pdf);
    randomize(field10, engine, pdf);
    randomize(field11, engine, pdf);
    randomize(field20, engine, pdf);

    FieldIdentifier F_x("fx", {{COL}, {dim0}}, m / s, "g");
    FieldIdentifier F_y("fy", {{LEV}, {dim2}}, m / s, "g");
    FieldIdentifier F_z("fz", {{CMP}, {dim1}}, m / s, "g");
    FieldIdentifier F_w("fyz", {{CMP, LEV}, {dim1, dim2}}, m / s, "g");

    Field field_x(F_x);
    Field field_y(F_y);
    Field field_z(F_z);
    Field field_w(F_w);

    // Test invalid inputs
    REQUIRE_THROWS(horiz_contraction<Real>(fieldsc, field_x,
                                           field00));  // x not allocated yet

    field_x.allocate_view();
    field_y.allocate_view();
    field_z.allocate_view();
    field_w.allocate_view();

    REQUIRE_THROWS(horiz_contraction<Real>(fieldsc, field_y,
                                           field_x));  // unmatching layout
    REQUIRE_THROWS(horiz_contraction<Real>(field_z, field11,
                                           field11));  // wrong weight layout

    Field result;

    // Ensure a scalar case works
    result = fieldsc.clone();
    horiz_contraction<Real>(result, field00, field00);
    result.sync_to_host();
    auto v = result.get_view<Real, Host>();
    // The numerator is the sum of the squares of the first dim0 integers
    // The denominator is the sum of the first dim0 integers squared
    Real wavg = sp(sum_n_sq(dim0)) / sp(sum_n(dim0) * sum_n(dim0));
    REQUIRE_THAT(v(), Catch::Matchers::WithinRel(wavg, tol));

    // Test higher-order cases
    result = field_z.clone();
    horiz_contraction<Real>(result, field10, field00);
    REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
            std::vector<FieldTag>({CMP}));
    REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim1);

    result = field_y.clone();
    horiz_contraction<Real>(result, field11, field00);
    REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
            std::vector<FieldTag>({LEV}));
    REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim2);

    result = field_w.clone();
    horiz_contraction<Real>(result, field20, field00);
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
      using RPDF  = std::uniform_real_distribution<Real>;
      auto engine = setup_random_test();
      RPDF pdf(0, 1);

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
      randomize(fieldsc, engine, pdf);
      randomize(field10, engine, pdf);
      randomize(field11, engine, pdf);
      randomize(field20, engine, pdf);

      FieldIdentifier F_x("fx", {{COL}, {dim0}}, m / s, "g");
      FieldIdentifier F_y("fy", {{CMP}, {dim1}}, m / s, "g");
      FieldIdentifier F_z("fz", {{COL, CMP}, {dim0, dim1}}, m / s, "g");

      Field field_x(F_x);
      Field field_y(F_y);
      Field field_z(F_z);

      // Test invalid inputs
      REQUIRE_THROWS(vert_contraction<Real>(fieldsc, field_x,
                                            field00));  // x not allocated yet

      field_x.allocate_view();
      field_y.allocate_view();
      field_z.allocate_view();

      REQUIRE_THROWS(vert_contraction<Real>(fieldsc, field_y,
                                            field_x));  // unmatching layout
      REQUIRE_THROWS(vert_contraction<Real>(field_z, field11,
                                            field11));  // wrong weight layout

      Field result;

      // Add test for invalid rank-2 weight field layout
      FieldIdentifier bad_w("bad_w", {{CMP, lev_tag}, {dim1, dim2}}, m / s,
                            "g");
      Field bad_weight(bad_w);
      bad_weight.allocate_view();
      REQUIRE_THROWS(vert_contraction<Real>(result, field20, bad_weight));

      // Add test for mismatched weight field dimensions
      FieldIdentifier wrong_size_w(
          "wrong_w", {{COL, lev_tag}, {dim0 + 1, dim2}}, m / s, "g");
      Field wrong_weight(wrong_size_w);
      wrong_weight.allocate_view();
      REQUIRE_THROWS(vert_contraction<Real>(result, field20, wrong_weight));

      // Ensure a scalar case works
      result = fieldsc.clone();
      vert_contraction<Real>(result, field00, field00);
      result.sync_to_host();
      auto v = result.get_view<Real, Host>();
      // The numerator is the sum of the squares of the first dim2 integers
      // (analytically known). The denominator is the sum of the first dim2
      // integers squared (analytically known)
      Real havg = sp(sum_n_sq(dim2)) / sp(sum_n(dim2) * sum_n(dim2));
      REQUIRE_THAT(v(), Catch::Matchers::WithinRel(havg, tol));

      // Test higher-order cases
      result = field_x.clone();
      vert_contraction<Real>(result, field10, field00);
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
      vert_contraction<Real>(result, field11, field00);
      REQUIRE(result.get_header().get_identifier().get_layout().tags() ==
              std::vector<FieldTag>({CMP}));
      REQUIRE(result.get_header().get_identifier().get_layout().dim(0) == dim1);

      result = field_z.clone();
      vert_contraction<Real>(result, field20, field00);
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
      vert_contraction<Real>(result, field20, field10);
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

  SECTION ("frobenius") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    // Summing the squares from n=a+1 to n=a+N gives
    // (a+1)^2+(a+2)^2+...
    // which ultimately gives
    // N*a^2 + 2a*(1+2+...+N)
    Real lsum = offset*offset*lsize + 2*offset*lsize*(lsize+1)/2 + lsize*(lsize+1)*(2*lsize+1) / 6.0;
    Real gsum = gsize*(gsize+1)*(2*gsize+1) / 6.0;

    REQUIRE(frobenius_norm<Real>(f1)==std::sqrt(lsum));
    REQUIRE(frobenius_norm<Real>(f1,&comm)==std::sqrt(gsum));
  }

  SECTION ("max") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmax = offset + lsize;
    Real gmax = gsize;

    REQUIRE(field_max<Real>(f1)==lmax);
    REQUIRE(field_max<Real>(f1,&comm)==gmax);
  }

  SECTION ("min") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmin = offset + 1;
    Real gmin = 1;

    REQUIRE(field_min<Real>(f1)==lmin);
    REQUIRE(field_min<Real>(f1,&comm)==gmin);
  }

  SECTION ("perturb") {
    using namespace ShortFieldTagsNames;
    using RPDF = std::uniform_real_distribution<Real>;
    using IPDF = std::uniform_int_distribution<int>;
    auto engine = setup_random_test ();

    const int ncols = 6;
    const int ncmps = 2;
    const int nlevs = IPDF(3,9)(engine); // between 3-9 levels

    // Create 1d, 2d, 3d fields with a level dimension, and set all to 1
    FieldIdentifier fid1 ("f_1d",   FieldLayout({LEV},           {nlevs}),               Units::nondimensional(), "");
    FieldIdentifier fid2a("f_2d_a", FieldLayout({CMP, LEV},      {ncmps, nlevs}),        Units::nondimensional(), "");
    FieldIdentifier fid2b("f_2d_b", FieldLayout({COL, LEV},      {ncols, nlevs}),        Units::nondimensional(), "");
    FieldIdentifier fid3 ("f_3d",   FieldLayout({COL, CMP, LEV}, {ncols, ncmps, nlevs}), Units::nondimensional(), "");
    Field f1(fid1), f2a(fid2a), f2b(fid2b), f3(fid3);
    f1.allocate_view(), f2a.allocate_view(), f2b.allocate_view(), f3.allocate_view();
    f1.deep_copy(1), f2a.deep_copy(1), f2b.deep_copy(1), f3.deep_copy(1);

    // We need GIDs for fields with COL component. This test is not over
    // multiple ranks, so just set as [0, ncols-1].
    Field gids(FieldIdentifier("gids", FieldLayout({COL}, {ncols}), Units::nondimensional(), "", DataType::IntType));
    gids.allocate_view();
    auto gids_data = gids.get_internal_view_data<int,Host>();
    std::iota(gids_data, gids_data+ncols, 0);
    gids.sync_to_dev();

    // Create masks s.t. only last 3 levels are perturbed. For variety,
    // 1d and 2d fields will use lambda mask and 3 field will use a view.
    auto mask_lambda = [&nlevs] (const int& i0) {
      return i0 >= nlevs-3;
    };
    kt_host::view_1d<bool> mask_view("mask_view", nlevs);
    Kokkos::deep_copy(mask_view, false);
    for (int ilev=0; ilev<nlevs; ++ilev) {
      if (ilev >= nlevs-3) mask_view(ilev) = true;
    }

    // Compute random perturbation between [2, 3]
    RPDF pdf(2, 3);
    int base_seed = 0;
    perturb(f1,  engine, pdf, base_seed, mask_lambda);
    perturb(f2a, engine, pdf, base_seed, mask_lambda);
    perturb(f2b, engine, pdf, base_seed, mask_lambda, gids);
    perturb(f3,  engine, pdf, base_seed, mask_view,   gids);

    // Sync to host for checks
    f1.sync_to_host(), f2a.sync_to_host(), f2b.sync_to_host(), f3.sync_to_host();
    const auto v1  = f1.get_strided_view <Real*,   Host>();
    const auto v2a = f2a.get_strided_view<Real**,  Host>();
    const auto v2b = f2b.get_strided_view<Real**,  Host>();
    const auto v3  = f3.get_strided_view <Real***, Host>();

    // Check that all field values are 1 for all but last 3 levels and between [2,3] otherwise.
    auto check_level = [&] (const int ilev, const Real val) {
      if (ilev < nlevs-3) REQUIRE(val == 1);
      else REQUIRE((2 <= val && val <= 3));
    };
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmps; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (icol==0 && icmp==0) check_level(ilev, v1(ilev));
          if (icol==0) check_level(ilev, v2a(icmp,ilev));
          if (icmp==0) check_level(ilev, v2b(icol,ilev));
          check_level(ilev, v3(icol,icmp,ilev));
    }}}

    // Check that using a different seed gives different values
    auto f1_alt = f1.clone(); f1_alt.deep_copy(1.0);
    auto f3_alt = f3.clone(); f3_alt.deep_copy(1.0);
    int base_seed_alt = 100;
    perturb(f1_alt, engine, pdf, base_seed_alt, mask_lambda);
    perturb(f3_alt, engine, pdf, base_seed_alt, mask_lambda, gids);
    f1_alt.sync_to_host(), f3_alt.sync_to_host();

    const auto v1_alt = f1_alt.get_strided_view<Real*,   Host>();
    const auto v3_alt = f3_alt.get_strided_view<Real***, Host>();

    auto check_diff = [&] (const int ilev, const Real val1, const Real val2) {
      if (ilev < nlevs-3) REQUIRE(val1==val2);
      else                REQUIRE(val1!=val2);
    };
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmps; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (icol==0 && icmp==0) check_diff(ilev, v1(ilev), v1_alt(ilev));
          check_diff(ilev, v3(icol,icmp,ilev), v3_alt(icol,icmp,ilev));
    }}}

    // Finally check that the original seed gives same result
    f1_alt.deep_copy(1.0), f3_alt.deep_copy(1.0);
    perturb(f1_alt, engine, pdf, base_seed, mask_lambda);
    perturb(f3_alt, engine, pdf, base_seed, mask_lambda, gids);
    REQUIRE(views_are_equal(f1, f1_alt));
    REQUIRE(views_are_equal(f3, f3_alt));
  }

  SECTION ("wrong_st") {
    using wrong_real =
      typename std::conditional<std::is_same<Real,double>::value,
                                float,double>::type;
    REQUIRE_THROWS(field_min<int>(f1));
    REQUIRE_THROWS(field_max<int>(f1));
    REQUIRE_THROWS(field_sum<wrong_real>(f1));
    REQUIRE_THROWS(frobenius_norm<wrong_real>(f1));
  }
}

TEST_CASE ("print_field_hyperslab") {
  using namespace scream;
  using namespace ekat::units;

  using namespace ShortFieldTagsNames;
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  auto engine = setup_random_test ();
  RPDF pdf(0,1);

  const int nel = 2;
  const int ncmp = 2;
  const int ngp  = 4;
  const int nlev = 12;
  const int iel  = IPDF(0,nel-1)(engine);
  const int icmp = IPDF(0,ncmp-1)(engine);
  const int igp  = IPDF(0,ngp-1)(engine);
  const int jgp  = IPDF(0,ngp-1)(engine);
  const int ilev = IPDF(0,nlev-1)(engine);

  // WARNING: make sure this value matches the one used in print_field_hyperslab
  constexpr int max_per_line = 5;

  // Create field (if available, use packs, to ensure we don't print garbage)
  std::vector<FieldTag> tags = {EL, CMP, GP, GP, LEV};
  std::vector<int>      dims = {nel,ncmp,ngp,ngp,nlev};
  int pack_size = std::min(SCREAM_PACK_SIZE,4);

  FieldIdentifier fid ("f", {tags,dims}, kg, "some_grid");
  Field f (fid);
  f.get_header().get_alloc_properties().request_allocation(pack_size);
  f.allocate_view();
  randomize (f,engine,pdf);

  auto v = f.get_strided_view<const Real*****,Host>();

  SECTION ("slice_0") {
    std::vector<FieldTag> loc_tags = {EL,CMP};
    std::vector<int>      loc_idxs = {iel,icmp};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
      for (int gp1=0; gp1<ngp; ++gp1) {
        for (int gp2=0; gp2<ngp; ++gp2) {
          expected <<  "  f(" << iel << "," << icmp << "," << gp1 << "," << gp2 << ",:)";
          for (int lev=0; lev<nlev; ++lev) {
            if (lev % max_per_line == 0) {
              expected << "\n    ";
            }
            expected << v(iel,icmp,gp1,gp2,lev) << ", ";
          }
          expected << "\n";
        }
      }

    REQUIRE (out.str()==expected.str());
  }
  SECTION ("slice_0234") {
    std::vector<FieldTag> loc_tags = {EL,GP,GP,LEV};
    std::vector<int>      loc_idxs = {iel,igp,jgp,ilev};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
      expected <<  "  f(" << iel << ",:," << igp << "," << jgp << "," << ilev << ")";
      for (int cmp=0; cmp<ncmp; ++cmp) {
        if (cmp % max_per_line == 0) {
          expected << "\n    ";
        }
        expected << v(iel,cmp,igp,jgp,ilev) << ", ";
      }
      expected << "\n";

    REQUIRE (out.str()==expected.str());
  }
  SECTION ("slice_01234") {
    std::vector<FieldTag> loc_tags = {EL,CMP,GP,GP,LEV};
    std::vector<int>      loc_idxs = {iel,icmp,igp,jgp,ilev};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
      expected <<  "  f(" << ekat::join(loc_idxs,",") << ")\n";
      expected << "    " << v(iel,icmp,igp,jgp,ilev) << ", \n";

    REQUIRE (out.str()==expected.str());
  }
}

TEST_CASE ("compute_mask") {
  using namespace scream;

  using namespace ShortFieldTagsNames;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);

  const int ncols = 3*comm.size();
  const int nlevs = 128;
  const auto units = ekat::units::Units::nondimensional();

  // Create field (if available, use packs, to ensure we don't print garbage)
  std::vector<FieldTag> tags3d = {COL, CMP, LEV};
  std::vector<FieldTag> tags2d = {COL, LEV};
  std::vector<int>      dims3d = {ncols,2,nlevs};
  std::vector<int>      dims2d = {ncols,nlevs};

  FieldIdentifier fid3d ("foo", {tags3d,dims3d}, units, "some_grid");
  FieldIdentifier fid3di ("foo", {tags3d,dims3d}, units, "some_grid", DataType::IntType);
  FieldIdentifier fid2d ("foo", {tags2d,dims2d}, units, "some_grid");

  SECTION ("exceptions") {
    // Test compute_mask exception handling
    Field f (fid3d);
    Field m1 (fid3d);

    REQUIRE_THROWS(compute_mask<Comparison::EQ>(f,1,m1)); // Field not allocated
    f.allocate_view();
    REQUIRE_THROWS(compute_mask<Comparison::EQ>(f,1,m1)); // Mask not allocated
    m1.allocate_view();

    Field m2 (fid2d);
    m2.allocate_view();
    REQUIRE_THROWS(compute_mask<Comparison::EQ>(f,1,m2)); // incompatible layouts
  }

  SECTION ("check") {
    Field x(fid3d), one(fid3di), zero(fid3di), m(fid3di);

    x.allocate_view();
    one.allocate_view();
    m.allocate_view();
    zero.allocate_view();

    one.deep_copy(1);
    zero.deep_copy(0);
    x.deep_copy(2);

    // x==1 is false
    m.deep_copy(-1);
    compute_mask<Comparison::EQ>(x,1,m);
    REQUIRE(views_are_equal(m,zero));

    // x!=1 is true
    m.deep_copy(-1);
    compute_mask<Comparison::NE>(x,1,m);
    REQUIRE(views_are_equal(m,one));

    // x==2 is true
    m.deep_copy(-1);
    compute_mask<Comparison::EQ>(x,2,m);
    REQUIRE(views_are_equal(m,one));

    // x>1 is true
    m.deep_copy(-1);
    compute_mask<Comparison::GT>(x,1,m);
    REQUIRE(views_are_equal(m,one));

    // x>2 is false
    m.deep_copy(-1);
    compute_mask<Comparison::GT>(x,2,m);
    REQUIRE(views_are_equal(m,zero));

    // x>=2 is true
    m.deep_copy(-1);
    compute_mask<Comparison::GE>(x,2,m);
    REQUIRE(views_are_equal(m,one));

    // x<3 is true
    m.deep_copy(-1);
    compute_mask<Comparison::LT>(x,3,m);
    REQUIRE(views_are_equal(m,one));

    // x<2 is flase
    m.deep_copy(-1);
    compute_mask<Comparison::LT>(x,2,m);
    REQUIRE(views_are_equal(m,zero));

    // x<=2 is true
    m.deep_copy(-1);
    compute_mask<Comparison::LE>(x,2,m);
    REQUIRE(views_are_equal(m,one));
  }
}

} // anonymous namespace
