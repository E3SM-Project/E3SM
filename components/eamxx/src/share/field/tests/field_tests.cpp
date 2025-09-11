#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/field/field_impl.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>

namespace {

TEST_CASE("field", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  using P4 = ekat::Pack<Real,4>;
  using P8 = ekat::Pack<Real,8>;
  using P16 = ekat::Pack<Real,16>;

  auto engine = setup_random_test ();
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.01,0.99);

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");

  // Check if we can extract a reshaped view
  SECTION ("reshape") {
    Field f1 (fid);

    // Should not be able to reshape before allocating
    REQUIRE_THROWS(f1.get_view<Real*>());

    f1.allocate_view();

    // Reshape should work with both dynamic and static dims
    auto v1 = f1.get_view<Real[3][24]>();
    auto v2 = f1.get_view<Real**>();

    REQUIRE(v1.size()==v2.size());

    // But if wrong static length is used, we should throw
    REQUIRE_THROWS(f1.get_view<Real[3][16]>());

    // Should not be able to reshape to this data type...
    REQUIRE_THROWS(f1.get_view<P16**>());
    // But this should work
    f1.get_view<P8**>();

    // Using packs (of allowable size) of different pack sizes
    // should lead to views with different extents.
    // Since there's no padding, their extent on last dimension
    // should be the phys dim divided by pack size.
    auto v3 = f1.get_view<P8**>();
    auto v4 = f1.get_view<P4**>();
    REQUIRE (v4.size() == 2*v3.size());
    REQUIRE (v4.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v3.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v4.extent_int(1) == fid.get_layout().dim(1) / P4::n);
    REQUIRE (v3.extent_int(1) == fid.get_layout().dim(1) / P8::n);

    // The memory spans should be identical
    REQUIRE (v3.impl_map().memory_span()==v4.impl_map().memory_span());

    // Trying to reshape into something that the allocation cannot accommodate should throw
    REQUIRE_THROWS (f1.get_view<P16***>());

    // Can't get non-const data type view from a read-only field
    REQUIRE_THROWS (f1.get_const().get_view<Real**>());
  }

  // Check copy constructor
  SECTION ("copy ctor") {
    Field f1 (fid);

    f1.allocate_view();
    f1.deep_copy(3.0);

    Field f2 = f1;
    REQUIRE(f2.get_header_ptr()==f1.get_header_ptr());
    REQUIRE(f2.get_internal_view_data<Real>()==f1.get_internal_view_data<Real>());
    REQUIRE(f2.is_allocated());
    REQUIRE(views_are_equal(f1,f2));
  }

  SECTION ("construct_from_view") {
    // Crate f1 with some padding, to stress test the feature
    Field f1 (fid);
    auto& fap1 = f1.get_header().get_alloc_properties();
    fap1.request_allocation(16);
    f1.allocate_view();
    f1.deep_copy(1.0);

    // Get f1 view, and wrap it in another field
    auto view = f1.get_view<Real**>();
    Field f2 (fid,view);

    // Check the two are the same
    REQUIRE (views_are_equal(f1,f2));

    // Modify one field, and check again
    randomize(f2,engine,pdf);
    REQUIRE (views_are_equal(f1,f2));
  }

  SECTION ("clone") {
    Field f1 (fid);
    auto& fap1 = f1.get_header().get_alloc_properties();

    fap1.request_allocation(16);
    f1.allocate_view();
    f1.deep_copy(3.0);

    Field f2 = f1.clone();
    auto& fap2 = f2.get_header().get_alloc_properties();
    REQUIRE(f2.is_allocated());
    REQUIRE(fap2.get_alloc_size()==fap1.get_alloc_size());
    REQUIRE(views_are_equal(f1,f2));

    // Changing f2 should leave f1 unchanged
    f2.deep_copy(0);
    REQUIRE (field_max<Real>(f2)==0.0);
    REQUIRE (field_min<Real>(f2)==0.0);
    REQUIRE (field_max<Real>(f1)==3.0);
    REQUIRE (field_min<Real>(f1)==3.0);
  }

  SECTION ("alias") {
    Field f1 (fid);
    f1.allocate_view();

    Field f2 = f1.alias("the_alias");

    REQUIRE(f2.is_allocated());
    REQUIRE(&f1.get_header().get_tracking()==&f2.get_header().get_tracking());
    REQUIRE(&f1.get_header().get_alloc_properties()==&f2.get_header().get_alloc_properties());
    REQUIRE(f1.get_header().get_identifier().get_layout()==f2.get_header().get_identifier().get_layout());
    REQUIRE(f1.get_internal_view_data<Real>()==f2.get_internal_view_data<Real>());

    // Identifiers are separate objects though
    REQUIRE(&f1.get_header().get_identifier()!=&f2.get_header().get_identifier());

    // Check extra data is also shared
    f1.get_header().set_extra_data("foo",1);
    REQUIRE (f2.get_header().has_extra_data("foo"));
  }

  SECTION ("is_aliasing") {
    Field f1 (fid), f2(fid);
    f1.allocate_view();
    f2.allocate_view();

    // Check aliasing
    REQUIRE (f1.is_aliasing(f1.get_const()));
    REQUIRE (f1.is_aliasing(f1));
    REQUIRE (f1.is_aliasing(f1.alias("foo")));

    // f1 and f2 have independent views, so they are not aliasing each other
    REQUIRE (!f1.is_aliasing(f2));
    REQUIRE (!f1.is_aliasing(f2));

    auto f1_0x = f1.subfield(0,0);
    auto f1_1x = f1.subfield(0,1);
    auto f1_x0 = f1.subfield(1,0);
    auto f1_x1 = f1.subfield(1,1);

    auto g1_0x = f1.subfield(0,0);
    auto g1_1x = f1.subfield(0,1);
    auto g1_x0 = f1.subfield(1,0);
    auto g1_x1 = f1.subfield(1,1);

    // Check we preserve parent info
    auto f1_0x_p = f1_0x.get_header().get_parent();
    REQUIRE (f1_0x.alias("foo").get_header().get_parent()==f1_0x_p);

    REQUIRE (f1_0x.is_aliasing(g1_0x));
    REQUIRE (f1_1x.is_aliasing(g1_1x));
    REQUIRE (f1_x0.is_aliasing(g1_x0));
    REQUIRE (f1_x1.is_aliasing(g1_x1));

    REQUIRE (not f1_0x.is_aliasing(f1_1x));
    REQUIRE (not f1_0x.is_aliasing(f1_x0));
    REQUIRE (not f1_0x.is_aliasing(f1_x1));

    REQUIRE (not f1_1x.is_aliasing(f1_x0));
    REQUIRE (not f1_1x.is_aliasing(f1_x1));

    REQUIRE (not f1_x0.is_aliasing(f1_x1));
  }

  SECTION ("deep_copy") {
    // rank-0
    std::vector<FieldTag> t0 = {};
    std::vector<int> d0 = {};
    FieldIdentifier fid0("scalar_0d",{t0,d0},m/s,"some_grid");
    Field f0(fid0);
    f0.allocate_view();
    f0.deep_copy(1.5);
    f0.sync_to_host();
    REQUIRE (reinterpret_cast<Real*>(f0.get_internal_view_data<Real,Host>())[0]==1.5);

    // rank-3
    std::vector<FieldTag> t1 = {COL,CMP,LEV};
    std::vector<int> d1 = {3,2,24};

    FieldIdentifier fid1("vec_3d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    f1.deep_copy(1.0);
    f1.sync_to_host();
    auto v = reinterpret_cast<Real*>(f1.get_internal_view_data<Real,Host>());
    for (int i=0; i<fid1.get_layout().size(); ++i) {
      REQUIRE (v[i]==1.0);
    }
  }

  SECTION ("host_view") {
    Field f(fid);

    // Views not yet allocated
    REQUIRE_THROWS(f.get_internal_view_data<Real>());
    REQUIRE_THROWS(f.get_internal_view_data<Real,Host>());
    REQUIRE_THROWS(f.sync_to_host());
    REQUIRE_THROWS(f.sync_to_dev());

    f.allocate_view();
    randomize(f,engine,pdf);

    // Get reshaped view on device, and manually create Host mirror
    auto v2d = f.get_view<Real**>();
    auto v2d_hm = Kokkos::create_mirror_view(v2d);
    Kokkos::deep_copy(v2d_hm,v2d);

    // Get reshaped view straight on Host
    auto v2dh = f.get_view<Real**,Host>();

    // The two should match
    for (int i=0; i<dims[0]; ++i) {
      for (int j=0; j<dims[1]; ++j) {
        REQUIRE (v2dh(i,j) == v2d_hm(i,j) );
      }
    }
  }

  SECTION ("rank0_field") {
    // Create 0d field
    FieldIdentifier fid0("f_0d", FieldLayout({},{}), Units::nondimensional(), "dummy_grid");
    Field f0(fid0);
    f0.allocate_view();

    // Create 1d field
    FieldIdentifier fid1("f_1d", FieldLayout({COL}, {5}), Units::nondimensional(), "dummy_grid");
    Field f1(fid1);
    f1.allocate_view();

    // Randomize 1d field
    randomize(f1,engine,pdf);

    auto v0 = f0.get_view<Real, Host>();
    auto v1 = f1.get_view<Real*, Host>();

    // Deep copy subfield of 1d field -> 0d field and check result
    for (size_t i=0; i<v1.extent(0); ++i) {
      f0.deep_copy<Host>(f1.subfield(0, i));
      REQUIRE(v0() == v1(i));
    }

    // Randomize 0d field
    randomize(f0,engine,pdf);

    // Deep copy 0d field -> subfield of 1d field and check result
    for (size_t i=0; i<v1.extent(0); ++i) {
      f1.subfield(0, i).deep_copy<Host>(f0);
      REQUIRE(v1(i) == v0());
    }
  }
}

} // anonymous namespace
