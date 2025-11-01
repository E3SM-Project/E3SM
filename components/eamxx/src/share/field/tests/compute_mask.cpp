#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

TEST_CASE ("compute_mask") {
  using namespace ShortFieldTagsNames;

  const int ncols = 10;
  const int nlevs = 128;
  const auto units = ekat::units::Units::nondimensional();

  // Create fields
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

    REQUIRE_THROWS(compute_mask(f,1,Comparison::EQ,m1)); // Field not allocated
    f.allocate_view();
    REQUIRE_THROWS(compute_mask(f,1,Comparison::EQ,m1)); // Mask not allocated
    m1.allocate_view();

    Field m2 (fid2d);
    m2.allocate_view();
    REQUIRE_THROWS(compute_mask(f,1,Comparison::EQ,m2)); // incompatible layouts
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
    compute_mask(x,1,Comparison::EQ,m);
    REQUIRE(views_are_equal(m,zero));

    // x!=1 is true
    m.deep_copy(-1);
    compute_mask(x,1,Comparison::NE,m);
    REQUIRE(views_are_equal(m,one));

    // x==2 is true
    m.deep_copy(-1);
    compute_mask(x,2,Comparison::EQ,m);
    REQUIRE(views_are_equal(m,one));

    // x>1 is true
    m.deep_copy(-1);
    compute_mask(x,1,Comparison::GT,m);
    REQUIRE(views_are_equal(m,one));

    // x>2 is false
    m.deep_copy(-1);
    compute_mask(x,2,Comparison::GT,m);
    REQUIRE(views_are_equal(m,zero));

    // x>=2 is true
    m.deep_copy(-1);
    compute_mask(x,2,Comparison::GE,m);
    REQUIRE(views_are_equal(m,one));

    // x<3 is true
    m.deep_copy(-1);
    compute_mask(x,3,Comparison::LT,m);
    REQUIRE(views_are_equal(m,one));

    // x<2 is flase
    m.deep_copy(-1);
    compute_mask(x,2,Comparison::LT,m);
    REQUIRE(views_are_equal(m,zero));

    // x<=2 is true
    m.deep_copy(-1);
    compute_mask(x,2,Comparison::LE,m);
    REQUIRE(views_are_equal(m,one));
  }
}

} // namespace scream
