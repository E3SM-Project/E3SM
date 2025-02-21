#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("field", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto engine = setup_random_test();
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.01, 0.99);

  // Subfields
  SECTION("subfield") {
    std::vector<FieldTag> t1 = {COL, CMP, CMP, LEV};
    std::vector<int> d1 = {3, 10, 2, 24};

    FieldIdentifier fid1("4d", {t1, d1}, m / s, "some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1, engine, pdf);

    const int idim = 1;
    const int ivar = 2;

    auto f2 = f1.subfield(idim, ivar);

    // Wrong rank for the subfield f2
    REQUIRE_THROWS(f2.get_view<Real****>());

    auto v4d_h = f1.get_view<Real****, Host>();
    auto v3d_h = f2.get_view<Real***, Host>();
    for (int i = 0; i < d1[0]; ++i)
      for (int j = 0; j < d1[2]; ++j)
        for (int k = 0; k < d1[3]; ++k) {
          REQUIRE(v4d_h(i, ivar, j, k) == v3d_h(i, j, k));
        }
  }

  SECTION("multi-sliced subfield") {
    SECTION("1D multi-slice") {
      // ==============
      /* Rank-1 view */
      // ==============
      std::vector<FieldTag> t1 = {COL};
      std::vector<int> d1 = {11};
      FieldIdentifier fid1("1d", {t1, d1}, m / s, "some_grid");

      Field f1(fid1);
      f1.allocate_view();
      randomize(f1, engine, pdf);

      const int idim = {0};
      const int sl_beg = {3};
      const int sl_end = {7};

      auto v1d_h = f1.get_view<Real*, Host>();
      auto sf = f1.subfield(idim, sl_beg, sl_end);
      REQUIRE(sf.get_header().get_alloc_properties().contiguous() == false);
      auto sv_h = sf.get_strided_view<Real*, Host>();
      REQUIRE(sv_h.extent_int(idim) == (sl_end - sl_beg));

      for (int i = sl_beg; i < sl_end; i++) {
        REQUIRE(v1d_h(i) == sv_h(i - sl_beg));
      }
    }

    SECTION("2D multi-slice") {
      std::vector<FieldTag> t2 = {COL, CMP};
      std::vector<int> d2 = {5, 10};
      FieldIdentifier fid2("2d", {t2, d2}, m / s, "some_grid");

      Field f2(fid2);
      f2.allocate_view();
      randomize(f2, engine, pdf);

      const int idim[2] = {0, 1};
      const int sl_beg[2] = {0, 3};
      const int sl_end[2] = {4, 6};

      auto v2d_h = f2.get_view<Real**, Host>();
      int i1, i2, j1, j2;

      for (int ens = 0; ens < 2; ens++) {
        auto sf = f2.subfield(idim[ens], sl_beg[ens], sl_end[ens]);
        auto sv_h = sf.get_strided_view<Real**, Host>();
        i1 = (ens == 0) ? sl_beg[0] : 0;
        i2 = (ens == 0) ? sl_end[0] : d2[0];
        j1 = (ens == 1) ? sl_beg[1] : 0;
        j2 = (ens == 1) ? sl_end[1] : d2[1];
        for (int i = i1; i < i2; i++) {
          for (int j = j1; j < j2; j++) {
            REQUIRE(v2d_h(i, j) == sv_h(i - i1, j - j1));
          }
        }
      }
    }
    SECTION("3D multi-slice") {
      // ==============
      /* Rank-3 view */
      // ==============
      std::vector<FieldTag> t3 = {COL, CMP, LEV};
      std::vector<int> d3 = {5, 10, 2};
      FieldIdentifier fid3("3d", {t3, d3}, m / s, "some_grid");

      Field f3(fid3);
      f3.allocate_view();
      randomize(f3, engine, pdf);

      const int idim[3] = {0, 1, 2};
      const int sl_beg[3] = {2, 3, 0};
      const int sl_end[3] = {4, 6, 1};

      auto v3d_h = f3.get_view<Real***, Host>();
      int i1, i2, j1, j2, k1, k2;

      for (int ens = 0; ens < 3; ens++) {
        auto sf = f3.subfield(idim[ens], sl_beg[ens], sl_end[ens]);
        auto sv_h = sf.get_strided_view<Real***, Host>();
        i1 = (ens == 0) ? sl_beg[0] : 0;
        i2 = (ens == 0) ? sl_end[0] : d3[0];
        j1 = (ens == 1) ? sl_beg[1] : 0;
        j2 = (ens == 1) ? sl_end[1] : d3[1];
        k1 = (ens == 2) ? sl_beg[2] : 0;
        k2 = (ens == 2) ? sl_end[2] : d3[2];
        for (int i = i1; i < i2; i++) {
          for (int j = j1; j < j2; j++) {
            for (int k = k1; k < k2; k++) {
              REQUIRE(v3d_h(i, j, k) == sv_h(i - i1, j - j1, k - k1));
            }
          }
        }
      }
      // ======================
      /*  get_components test */
      // ======================
      auto cmp3 = f3.get_components(sl_beg[1], sl_end[1]);

      auto svc_h = cmp3.get_strided_view<Real***, Host>();
      for (int i = 0; i < d3[0]; i++) {
        for (int j = sl_beg[1]; j < sl_end[1]; j++) {
          for (int k = 0; k < d3[2]; k++) {
            REQUIRE(v3d_h(i, j, k) == svc_h(i, j - sl_beg[1], k));
          }
        }
      }
    }
    SECTION("4D multi-slice") {
      // ==============
      /* Rank-4 view */
      // ==============
      std::vector<FieldTag> t4 = {COL, CMP, CMP, LEV};
      std::vector<int> d4 = {5, 10, 2, 23};
      FieldIdentifier fid4("4d", {t4, d4}, m / s, "some_grid");

      Field f4(fid4);
      f4.allocate_view();
      randomize(f4, engine, pdf);

      const int idim[4] = {0, 1, 2, 3};
      const int sl_beg[4] = {2, 3, 0, 9};
      const int sl_end[4] = {4, 6, 1, 15};

      auto v4d_h = f4.get_view<Real****, Host>();
      int i1, i2, j1, j2, k1, k2, l1, l2;

      for (int ens = 0; ens < 4; ens++) {
        auto sf = f4.subfield(idim[ens], sl_beg[ens], sl_end[ens]);
        auto sv_h = sf.get_strided_view<Real****, Host>();
        i1 = (ens == 0) ? sl_beg[0] : 0;
        i2 = (ens == 0) ? sl_end[0] : d4[0];
        j1 = (ens == 1) ? sl_beg[1] : 0;
        j2 = (ens == 1) ? sl_end[1] : d4[1];
        k1 = (ens == 2) ? sl_beg[2] : 0;
        k2 = (ens == 2) ? sl_end[2] : d4[2];
        l1 = (ens == 3) ? sl_beg[3] : 0;
        l2 = (ens == 3) ? sl_end[3] : d4[3];
        for (int i = i1; i < i2; i++) {
          for (int j = j1; j < j2; j++) {
            for (int k = k1; k < k2; k++) {
              for (int l = l1; l < l2; l++) {
                REQUIRE(v4d_h(i, j, k, l) ==
                        sv_h(i - i1, j - j1, k - k1, l - l1));
              }
            }
          }
        }
      }
    }
    SECTION("5D multi-slice") {
      // ==============
      /* Rank-5 view */
      // ==============
      std::vector<FieldTag> t5 = {EL, CMP, GP, GP, LEV};
      std::vector<int> d5 = {5, 10, 4, 2, 23};
      FieldIdentifier fid5("5d", {t5, d5}, m / s, "some_grid");

      Field f5(fid5);
      f5.allocate_view();
      randomize(f5, engine, pdf);

      const int idim[5] = {0, 1, 2, 3, 4};
      const int sl_beg[5] = {2, 3, 1, 0, 9};
      const int sl_end[5] = {4, 6, 3, 1, 15};

      auto v5d_h = f5.get_view<Real*****, Host>();
      int i1, i2, j1, j2, k1, k2, l1, l2, m1, m2;

      for (int ens = 0; ens < 5; ens++) {
        auto sf = f5.subfield(idim[ens], sl_beg[ens], sl_end[ens]);
        auto sv_h = sf.get_strided_view<Real*****, Host>();
        i1 = (ens == 0) ? sl_beg[0] : 0;
        i2 = (ens == 0) ? sl_end[0] : d5[0];
        j1 = (ens == 1) ? sl_beg[1] : 0;
        j2 = (ens == 1) ? sl_end[1] : d5[1];
        k1 = (ens == 2) ? sl_beg[2] : 0;
        k2 = (ens == 2) ? sl_end[2] : d5[2];
        l1 = (ens == 3) ? sl_beg[3] : 0;
        l2 = (ens == 3) ? sl_end[3] : d5[3];
        m1 = (ens == 4) ? sl_beg[4] : 0;
        m2 = (ens == 4) ? sl_end[4] : d5[4];
        for (int i = i1; i < i2; i++) {
          for (int j = j1; j < j2; j++) {
            for (int k = k1; k < k2; k++) {
              for (int l = l1; l < l2; l++) {
                for (int m = m1; m < m2; m++) {
                  REQUIRE(v5d_h(i, j, k, l, m) ==
                          sv_h(i - i1, j - j1, k - k1, l - l1, m - m1));
                }
              }
            }
          }
        }
      }
    }
    SECTION("6D multi-slice") {
      // ==============
      /* Rank-6 view */
      // ==============
      std::vector<FieldTag> t6 = {EL, TL, CMP, GP, GP, LEV};
      std::vector<int> d6 = {5, 10, 4, 2, 9, 23};
      FieldIdentifier fid6("6d", {t6, d6}, m / s, "some_grid");

      Field f6(fid6);
      f6.allocate_view();
      randomize(f6, engine, pdf);

      const int idim[6] = {0, 1, 2, 3, 4, 5};
      const int sl_beg[6] = {2, 3, 1, 0, 5, 9};
      const int sl_end[6] = {4, 6, 3, 1, 8, 15};

      auto v6d_h = f6.get_view<Real******, Host>();
      int i1, i2, j1, j2, k1, k2, l1, l2, m1, m2, n1, n2;

      for (int ens = 0; ens < 6; ens++) {
        auto sf = f6.subfield(idim[ens], sl_beg[ens], sl_end[ens]);
        REQUIRE(sf.get_header().get_alloc_properties().contiguous() == false);
        auto sv_h = sf.get_strided_view<Real******, Host>();
        i1 = (ens == 0) ? sl_beg[0] : 0;
        i2 = (ens == 0) ? sl_end[0] : d6[0];
        j1 = (ens == 1) ? sl_beg[1] : 0;
        j2 = (ens == 1) ? sl_end[1] : d6[1];
        k1 = (ens == 2) ? sl_beg[2] : 0;
        k2 = (ens == 2) ? sl_end[2] : d6[2];
        l1 = (ens == 3) ? sl_beg[3] : 0;
        l2 = (ens == 3) ? sl_end[3] : d6[3];
        m1 = (ens == 4) ? sl_beg[4] : 0;
        m2 = (ens == 4) ? sl_end[4] : d6[4];
        n1 = (ens == 5) ? sl_beg[5] : 0;
        n2 = (ens == 5) ? sl_end[5] : d6[5];
        for (int i = i1; i < i2; i++) {
          for (int j = j1; j < j2; j++) {
            for (int k = k1; k < k2; k++) {
              for (int l = l1; l < l2; l++) {
                for (int m = m1; m < m2; m++) {
                  for (int n = n1; n < n2; n++) {
                    REQUIRE(v6d_h(i, j, k, l, m, n) == sv_h(i - i1, j - j1,
                                                            k - k1, l - l1,
                                                            m - m1, n - n1));
                    REQUIRE((sv_h.extent_int(0) == (i2 - i1) &&
                             sv_h.extent_int(1) == (j2 - j1) &&
                             sv_h.extent_int(2) == (k2 - k1) &&
                             sv_h.extent_int(3) == (l2 - l1) &&
                             sv_h.extent_int(4) == (m2 - m1) &&
                             sv_h.extent_int(5) == (n2 - n1)));
                  }
                }
              }
            }
          }
        }
      }
    }
    SECTION("Subfield deep_copy") {
      // ==============
      /* Rank-3 view */
      // ==============
      std::vector<FieldTag> t3 = {COL, CMP, LEV};
      std::vector<int> d3 = {5, 10, 2};
      FieldIdentifier fid3("3d", {t3, d3}, m / s, "some_grid");

      Field f3a(fid3);
      f3a.allocate_view();
      randomize(f3a, engine, pdf);

      Field f3b(fid3);
      f3b.allocate_view();
      randomize(f3b, engine, pdf);

      const int idim = 1;
      const int sl_beg = 3;
      const int sl_end = 6;

      auto sfa = f3a.subfield(idim, sl_beg, sl_end);
      auto sv_a = sfa.get_strided_view<Real***, Host>();

      auto sfb = f3b.subfield(idim, sl_beg, sl_end);
      auto sv_b = sfb.get_strided_view<Real***, Host>();

      sfb.deep_copy<Host>(sfa);

      for (int i = 0; i < d3[0]; i++) {
        for (int j = sl_beg; j < sl_end; j++) {
          for (int k = 0; k < d3[2]; k++) {
            REQUIRE(sv_a(i, j - sl_beg, k) == sv_b(i, j - sl_beg, k));
          }
        }
      }
    }
  }
  // Dynamic Subfields
  SECTION("dynamic_subfield") {
    const int vec_dim = 10;
    std::vector<FieldTag> t1 = {COL, CMP, CMP, LEV};
    std::vector<int> d1 = {3, vec_dim, 2, 24};

    FieldIdentifier fid1("4d", {t1, d1}, m / s, "some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1, engine, pdf);

    const int idim = 1;
    const int ivar = 0;

    auto f2 = f1.subfield(idim, ivar, /* dynamic = */ true);

    // Cannot reset subview idx of non-subfield fields
    REQUIRE_THROWS(f1.get_header().get_alloc_properties().reset_subview_idx(0));

    // subview idx out of bounds
    auto& f2_ap = f2.get_header().get_alloc_properties();
    REQUIRE_THROWS(f2_ap.reset_subview_idx(-1));
    REQUIRE_THROWS(f2_ap.reset_subview_idx(vec_dim));

    // Fill f1 with random numbers, and verify corresponding subviews get same
    // values
    randomize(f1, engine, pdf);

    for (int ivar_dyn = 0; ivar_dyn < vec_dim; ++ivar_dyn) {
      // Reset slice idx
      f2_ap.reset_subview_idx(ivar_dyn);
      REQUIRE(f2_ap.get_subview_info().slice_idx == ivar_dyn);

      auto v4d_h = f1.get_view<Real****, Host>();
      auto v3d_h = f2.get_view<Real***, Host>();
      for (int i = 0; i < d1[0]; ++i)
        for (int j = 0; j < d1[2]; ++j)
          for (int k = 0; k < d1[3]; ++k) {
            REQUIRE(v4d_h(i, ivar_dyn, j, k) == v3d_h(i, j, k));
          }
    }
  }

  SECTION("vector_component") {
    std::vector<FieldTag> tags_2 = {COL, CMP, LEV};
    std::vector<int> dims_2 = {3, 2, 24};

    FieldIdentifier fid_2("vec_3d", {tags_2, dims_2}, m / s, "some_grid");

    Field f_vec(fid_2);
    f_vec.allocate_view();

    auto f0 = f_vec.get_component(0);
    auto f1 = f_vec.get_component(1);

    // No 3rd component
    REQUIRE_THROWS(f_vec.get_component(2));

    // f0 is scalar, no vector dimension
    REQUIRE_THROWS(f0.get_component(0));

    f0.deep_copy(1.0);
    f1.deep_copy(2.0);

    f_vec.sync_to_host();

    auto v = f_vec.get_view<Real***, Host>();
    for (int col = 0; col < 3; ++col) {
      for (int lev = 0; lev < 24; ++lev) {
        REQUIRE(v(col, 0, lev) == 1.0);
        REQUIRE(v(col, 1, lev) == 2.0);
      }
    }
  }
}
} // anonymous namespace
