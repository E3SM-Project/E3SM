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
    Field f2 = f_int.clone(CloneFlags::CopyData);

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

  SECTION ("scale") {
    SECTION ("real") {
      Field f1 = f_real.clone();
      Field f2 = f_real.clone(CloneFlags::CopyData);

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

  SECTION ("broadcast") {
    // f_real/f_int have layout (COL,CMP,LEV). Use rhs with fewer dims, broadcasted by the caller
    FieldLayout xl_cl ({COL,LEV},{ncol,nlev});
    FieldLayout xl_c  ({COL},{ncol});
    FieldLayout xl_l  ({LEV},{nlev});

    // Get x(i,j,k) for the broadcasted entries (using the original layout)
    auto xval = [&](const Field& x, const int i, const int k) {
      x.sync_to_host();
      const auto& l = x.get_header().get_identifier().get_layout();
      if (l.rank()==2) return x.get_strided_view<const Real**,Host>()(i,k);
      if (l.tags()[0]==COL) return x.get_strided_view<const Real*,Host>()(i);
      return x.get_strided_view<const Real*,Host>()(k);
    };

    for (const auto& xl : {xl_cl,xl_c,xl_l}) {
      Field x (FieldIdentifier("x",xl,kg,"some_grid"),true);
      randomize_uniform (x,seed++,1.0,2.0);
      auto xb = x.broadcast(f_real);

      f_real.sync_to_host();
      auto y0 = f_real.get_strided_view<const Real***,Host>();

      auto check = [&](const Field& y, auto&& expected) {
        y.sync_to_host();
        auto yv = y.get_strided_view<const Real***,Host>();
        for (int i=0; i<ncol; ++i)
          for (int j=0; j<ncmp; ++j)
            for (int k=0; k<nlev; ++k)
              REQUIRE_THAT (yv(i,j,k), Catch::Matchers::WithinRel(expected(i,j,k),1e-12));
      };

      { // scale
        Field y = f_real.clone(CloneFlags::CopyData);
        y.scale(xb);
        check (y,[&](int i,int j,int k){ return y0(i,j,k)*xval(x,i,k); });
      }
      { // scale_inv
        Field y = f_real.clone(CloneFlags::CopyData);
        y.scale_inv(xb);
        check (y,[&](int i,int j,int k){ return y0(i,j,k)/xval(x,i,k); });
      }
      { // update
        Field y = f_real.clone(CloneFlags::CopyData);
        y.update(xb,2.0,3.0);
        check (y,[&](int i,int j,int k){ return 2*xval(x,i,k)+3*y0(i,j,k); });
      }
      { // masked
        Field mask (fid_i.clone("mask"),true);
        mask.deep_copy(1);
        mask.subfield(COL,0).deep_copy(0);

        Field y = f_real.clone(CloneFlags::CopyData);
        y.scale(xb,mask);
        check (y,[&](int i,int j,int k){ return i==0 ? y0(i,j,k) : y0(i,j,k)*xval(x,i,k); });
      }
      { // src_mask
        // The broadcast of a masked field stores the broadcasted mask
        // Mask out the first entry of the first dim of x
        x.create_valid_mask(Field::MaskInit::Valid);
        x.get_valid_mask().subfield(0,0).deep_copy(0);
        auto xbm = x.broadcast(f_real);
        const bool first_is_col = xl.tags()[0]==COL;

        Field y = f_real.clone(CloneFlags::CopyData);
        y.scale(xbm,xbm.get_valid_mask());
        check (y,[&](int i,int j,int k){
          const bool masked_out = first_is_col ? i==0 : k==0;
          return masked_out ? y0(i,j,k) : y0(i,j,k)*xval(x,i,k);
        });
      }
    }

    // Int fields
    Field xi (FieldIdentifier("xi",xl_cl,kg,"some_grid",DataType::IntType),true);
    randomize_uniform (xi,seed++,1,5);
    Field yi = f_int.clone(CloneFlags::CopyData);
    yi.scale(xi.broadcast(f_int));
    yi.sync_to_host(); f_int.sync_to_host(); xi.sync_to_host();
    auto yiv = yi.get_strided_view<const int***,Host>();
    auto fiv = f_int.get_strided_view<const int***,Host>();
    auto xiv = xi.get_strided_view<const int**,Host>();
    for (int i=0; i<ncol; ++i)
      for (int j=0; j<ncmp; ++j)
        for (int k=0; k<nlev; ++k)
          REQUIRE (yiv(i,j,k)==fiv(i,j,k)*xiv(i,k));

    // Non-broadcasted rhs with different layout is still an error, and a broadcast cannot be the lhs
    Field x (FieldIdentifier("x",xl_cl,kg,"some_grid"),true);
    Field y = f_real.clone();
    REQUIRE_THROWS (y.scale(x));
    REQUIRE_THROWS (x.broadcast(f_real).scale(1.0));
  }

  SECTION ("max-min") {
    SECTION ("real") {
      Field one = f_real.clone();
      Field two = f_real.clone();
      one.deep_copy(1.0);
      two.deep_copy(2.0);

      Field f1 = one.clone(CloneFlags::CopyData);
      Field f2 = two.clone(CloneFlags::CopyData);
      f1.max(f2);
      REQUIRE (views_are_equal(f1, f2));

      Field f3 = one.clone(CloneFlags::CopyData);
      Field f4 = two.clone(CloneFlags::CopyData);
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

      Field f1 = one.clone(CloneFlags::CopyData);
      Field f2 = two.clone(CloneFlags::CopyData);
      f1.max(f2);
      REQUIRE (views_are_equal(f1, f2));

      Field f3 = one.clone(CloneFlags::CopyData);
      Field f4 = two.clone(CloneFlags::CopyData);
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
      Field f1 = f_real.clone(CloneFlags::CopyData);
      Field f2 = f_real.clone(CloneFlags::CopyData);
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
      Field f2 = f_real.clone(CloneFlags::CopyData);
      Field f3 = f_real.clone(CloneFlags::CopyData);

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
      REQUIRE (views_are_equal(f2,one));

      // Check handling of additive scalar
      f3.deep_copy(-2);
      f2.deep_copy(2);
      f3.update(f2,0,0,1);
      REQUIRE (views_are_equal(f3,one));
    }

    SECTION ("int") {
      Field f2 = f_int.clone(CloneFlags::CopyData);
      Field f3 = f_int.clone(CloneFlags::CopyData);

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

      // Check handling of additive scalar
      f3.deep_copy(-2);
      f2.deep_copy(2);
      f3.update(f2,0,0,1);
      REQUIRE (views_are_equal(f3,one));
    }
  }
}

} // anonymous namespace
