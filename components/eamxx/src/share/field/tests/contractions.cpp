#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

#include <limits>

namespace scream {

TEST_CASE("field_contractions") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // Accumulations in the Kokkos threaded reductions may be done in a
  // different order than the manual ones below, so we can only test
  // correctness up to a tolerance
  auto tol = std::numeric_limits<Real>::epsilon() * 100;

  auto seed = get_random_test_seed();

  SECTION("horiz_contraction") {

    int ncol = 129;
    int ncmp = 4;
    int nlev = 17;

    // Set a weight field with values [1,2,...,ncols]
    FieldIdentifier wfid("f", {{COL}, {ncol}}, none, "g");
    Field w(wfid,true);
    randomize_uniform(w,seed++);

    // Create (random) sample fields
    // NOTE: suffix _x stands for "horizontal", _z for "vertical", _xz for "horiz-and-vert"
    FieldIdentifier scl_fid("s", {{}, {}}, none, "g");  // scalar
    FieldIdentifier vec_x_fid("v_x", {{COL, CMP}, {ncol, ncmp}}, none, "g");
    FieldIdentifier scl_xz_fid("s_xz", {{COL, LEV}, {ncol, nlev}}, none, "g");
    FieldIdentifier vec_xz_fid("v_xz", {{COL, CMP, LEV}, {ncol, ncmp, nlev}}, none, "g");
    Field scl(scl_fid,true);
    Field vec_x(vec_x_fid,true);
    Field vec_xz(vec_xz_fid,true);
    randomize_uniform(scl, seed++);
    randomize_uniform(vec_x, seed++);
    randomize_uniform(vec_xz, seed++);

    FieldIdentifier scl_x_fid("s_x", {{COL}, {ncol}}, none, "g");
    FieldIdentifier scl_z_fid("s_z", {{LEV}, {nlev}}, none, "g");
    FieldIdentifier vec_fid  ("v", {{CMP}, {ncmp}}, none, "g");
    FieldIdentifier vec_z_fid("v_z", {{CMP, LEV}, {ncmp, nlev}}, none, "g");

    Field scl_x(scl_x_fid);
    Field vec  (vec_fid);
    Field vec_z(vec_z_fid);

    // Test invalid inputs
    REQUIRE_THROWS(horiz_contraction(scl, scl_x, scl_x.clone()));  // input not allocated yet
    REQUIRE_THROWS(horiz_contraction(vec, vec_x, scl_x.clone()));  // output not allocated yet
    REQUIRE_THROWS(horiz_contraction(scl, scl_x.clone(), scl_x));  // weight not allocated yet

    scl_x.allocate_view();
    vec.allocate_view();
    vec_z.allocate_view();

    REQUIRE_THROWS(horiz_contraction(scl, vec_x, scl_x)); // incompatible input-output layout
    REQUIRE_THROWS(horiz_contraction(scl, scl_x, vec_x)); // incompatible weight layout

    SECTION ("scl_x") {
      auto f_in = scl_x.clone();
      auto f_out = scl.clone();

      horiz_contraction(f_out, f_in, w);

      f_out.sync_to_host();
      auto v_in  = f_in.get_view<const Real *, Host>();
      auto v_out = f_out.get_view<Real, Host>();
      auto v_w = w.get_view<const Real*, Host>();
      Real manual = 0;
      for(int i = 0; i < ncol; ++i) {
        manual += v_w(i) * v_in(i);
      }
      REQUIRE_THAT(v_out(), Catch::Matchers::WithinRel(manual, tol));
    }

    SECTION ("scl_x_masked") {
      auto f_in = scl_x.clone();
      auto f_out = scl.clone();
      auto mask = f_in.create_valid_mask("mask",Field::MaskInit::Valid);
      mask.sync_to_host();
      auto v_mask = mask.get_view<int *, Host>();
      v_mask(ncol - 1) = 0;
      mask.sync_to_dev();
      
      horiz_contraction(f_out, f_in, w);
      f_out.sync_to_host();
      auto v_in  = f_in.get_view<const Real *, Host>();
      auto v_out = f_out.get_view<Real, Host>();
      auto v_w = w.get_view<const Real*, Host>();
      Real manual = 0;
      for(int i = 0; i < ncol; ++i) {
        if (v_mask(i))
          manual += v_w(i) * v_in(i);
      }
      REQUIRE_THAT(v_out(), Catch::Matchers::WithinRel(manual, tol));

      // Repeat but with ALL masked values (all mask=0 => sum = 0)
      mask.deep_copy(0);
      horiz_contraction(f_out, f_in, w);
      f_out.sync_to_host();
      REQUIRE_THAT(v_out(), Catch::Matchers::WithinAbs(Real(0), tol));
    }

    SECTION ("vec_xz") {
      auto f_in  = vec_xz.clone();
      auto f_out = vec_z.clone();
      auto mask = f_in.create_valid_mask();
      randomize_uniform(mask,seed++);

      horiz_contraction(f_out,f_in,w);

      f_out.sync_to_host();
      auto v_in = f_in.get_view<Real ***, Host>();
      auto v_out = f_out.get_view<Real **, Host>();
      auto v_mask = mask.get_view<const int ***, Host>();
      auto v_w = w.get_view<const Real*, Host>();

      for(int k = 0; k < nlev; ++k) {
        for(int j = 0; j < ncmp; ++j) {
          Real manual = 0;
          for(int i = 0; i < ncol; ++i) {
            if (v_mask(i,j,k))
              manual += v_w(i) * v_in(i, j, k);
          }
          REQUIRE_THAT(v_out(j, k), Catch::Matchers::WithinRel(manual, tol));
        }
      }
    }
  }

  SECTION("vert_contraction") {
    int ncol = 129;
    int ncmp = 4;
    int nlev = 17;

    // Test inputs with both LEV and ILEV vert dim tags
    for(FieldTag lev_tag : {LEV,ILEV}) {
      // Create (random) sample fields
      // NOTE: suffix _x stands for "horizontal", _z for "vertical", _xz for "horiz-and-vert"
      FieldIdentifier scl_fid   ("s",    {{}, {}}, ekat::units::none, "g");  // scalar
      FieldIdentifier scl_xz_fid("s_xz", {{COL, lev_tag}, {ncol, nlev}}, ekat::units::none, "g");
      FieldIdentifier vec_z_fid ("v_z",  {{CMP, lev_tag}, {ncmp, nlev}}, ekat::units::none, "g");
      FieldIdentifier vec_xz_fid("v_xz", {{COL, CMP, lev_tag}, {ncol, ncmp, nlev}}, ekat::units::none, "g");

      Field scl(scl_fid,true);
      Field scl_xz(scl_xz_fid,true);
      Field vec_z (vec_z_fid,true);
      Field vec_xz(vec_xz_fid,true);
      randomize_uniform(scl, seed++);
      randomize_uniform(scl_xz, seed++);
      randomize_uniform(vec_z, seed++);
      randomize_uniform(vec_xz, seed++);

      // Output/utility fields
      FieldIdentifier scl_x_fid("s_x", {{COL}, {ncol}}, ekat::units::none, "g");
      FieldIdentifier vec_x_fid("v_x", {{COL,CMP}, {ncol,ncmp}}, ekat::units::none, "g");
      FieldIdentifier vec_fid  ("v", {{CMP}, {ncmp}}, ekat::units::none, "g");
      FieldIdentifier scl_z_fid ("s_z",  {{lev_tag}, {nlev}}, ekat::units::none, "g");

      Field scl_x(scl_x_fid,true);
      Field vec_x(vec_x_fid,true);
      Field scl_z(scl_z_fid,true);
      Field vec(vec_fid,true);

      // Test invalid inputs
      REQUIRE_THROWS(vert_contraction(scl_x, Field(vec_z_fid), scl_z));  // input not allocated yet
      REQUIRE_THROWS(vert_contraction(Field(scl_x_fid), vec_z, scl_z));  // output not allocated yet
      REQUIRE_THROWS(vert_contraction(scl_x, vec_z, Field(scl_z_fid)));  // weight not allocated yet

      REQUIRE_THROWS(vert_contraction(scl, vec_z, scl_z));  // incompatible input-output layout
      REQUIRE_THROWS(vert_contraction(scl, scl_z, vec_xz)); // incompatible weight layout

      SECTION ("scl_z") {
        auto f_in = scl_z.clone();
        auto f_out = scl.clone();
        auto w = scl_z.clone();
        randomize_uniform (w,seed++);

        vert_contraction(f_out,f_in,w);
        f_out.sync_to_host();

        auto v_in  = f_in.get_view<const Real*,Host>(); 
        auto v_out = f_out.get_view<const Real,Host>(); 
        auto v_w   = w.get_view<const Real*,Host>(); 

        Real manual = 0;
        for (int lev=0; lev<nlev; ++lev) {
          manual += v_w(lev)*v_in(lev);
        }
        REQUIRE_THAT(v_out(), Catch::Matchers::WithinRel(manual, tol));
      }

      SECTION ("scl_z_masked") {
        auto f_in = scl_z.clone();
        auto f_out = scl.clone();
        auto w = scl_z.clone();
        randomize_uniform (w,seed++);
        auto mask = f_in.create_valid_mask();
        randomize_uniform (mask,seed++);

        vert_contraction(f_out,f_in,w);
        f_out.sync_to_host();

        auto v_in   = f_in.get_view<const Real*,Host>(); 
        auto v_out  = f_out.get_view<const Real,Host>(); 
        auto v_w    = w.get_view<const Real*,Host>(); 
        auto v_mask = mask.get_view<const int*,Host>(); 

        Real manual = 0;
        for (int lev=0; lev<nlev; ++lev) {
          if (v_mask(lev)) {
            manual += v_w(lev)*v_in(lev);
          }
        }
        REQUIRE_THAT(v_out(), Catch::Matchers::WithinRel(manual, tol));

        mask.deep_copy(0);
        vert_contraction(f_out,f_in,w);
        f_out.sync_to_host();
        REQUIRE_THAT(v_out(), Catch::Matchers::WithinRel(0, tol));
      }

      SECTION ("scl_xz_w_z") {
        auto f_in = scl_xz.clone();
        auto f_out = scl_x.clone();
        auto w = scl_z.clone();
        randomize_uniform (w,seed++);
        auto mask = f_in.create_valid_mask();
        randomize_uniform (mask,seed++);

        vert_contraction(f_out,f_in,w);
        f_out.sync_to_host();

        auto v_in   = f_in.get_view<const Real**,Host>(); 
        auto v_out  = f_out.get_view<const Real*,Host>(); 
        auto v_w    = w.get_view<const Real*,Host>(); 
        auto v_mask = mask.get_view<const int**,Host>(); 

        for (int col=0; col<ncol; ++col) {
          Real manual = 0;
          for (int lev=0; lev<nlev; ++lev) {
            if (v_mask(col,lev)) {
              manual += v_w(lev)*v_in(col,lev);
            }
          }
          REQUIRE_THAT(v_out(col), Catch::Matchers::WithinRel(manual, tol));
        }
      }

      SECTION ("vec_xz_w_xz") {
        auto f_in = vec_xz.clone();
        auto f_out = vec_x.clone();
        auto w = scl_xz.clone();
        randomize_uniform (w,seed++);
        auto mask = f_in.create_valid_mask();
        randomize_uniform (mask,seed++);

        vert_contraction(f_out,f_in,w);
        f_out.sync_to_host();

        auto v_in   = f_in.get_view<const Real***,Host>(); 
        auto v_out  = f_out.get_view<const Real**,Host>(); 
        auto v_w    = w.get_view<const Real**,Host>(); 
        auto v_mask = mask.get_view<const int***,Host>(); 

        for (int col=0; col<ncol; ++col) {
          for (int cmp=0; cmp<ncmp; ++cmp) {
            Real manual = 0;
            for (int lev=0; lev<nlev; ++lev) {
              if (v_mask(col,cmp,lev)) {
                manual += v_w(col,lev)*v_in(col,cmp,lev);
              }
            }
            REQUIRE_THAT(v_out(col,cmp), Catch::Matchers::WithinRel(manual, tol));
          }
        }
      }
    }
  }
}

} // namespace scream
