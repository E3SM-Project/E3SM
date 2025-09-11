#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/field/field_impl.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_utils.hpp>
#include <ekat_test_utils.hpp>
#include <ekat_subview_utils.hpp>

namespace {

TEST_CASE("field_layout", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;

  using TVec = std::vector<FieldTag>;
  using IVec = std::vector<int>;

  FieldLayout fl1 ({COL},{1});
  FieldLayout fl2 ({COL,CMP},{1,1});
  FieldLayout fl3 ({COL,CMP,CMP},{1,3,4});
  FieldLayout fl4 ({COL,LEV},{1,1});
  FieldLayout fl5 ({COL,CMP,LEV},{1,2,3});
  FieldLayout fl6 ({COL,CMP,CMP,ILEV},{1,5,6,1});
  FieldLayout fl7 ({LEV,CMP,COL},{3,2,1});

  REQUIRE (fl1.type()==LayoutType::Scalar2D);
  REQUIRE (fl2.type()==LayoutType::Vector2D);
  REQUIRE (fl3.type()==LayoutType::Tensor2D);
  REQUIRE (fl4.type()==LayoutType::Scalar3D);
  REQUIRE (fl5.type()==LayoutType::Vector3D);
  REQUIRE (fl6.type()==LayoutType::Tensor3D);

  REQUIRE (not fl1.is_vector_layout());
  REQUIRE (    fl2.is_vector_layout());
  REQUIRE (not fl3.is_vector_layout());
  REQUIRE (not fl4.is_vector_layout());
  REQUIRE (    fl5.is_vector_layout());
  REQUIRE (not fl6.is_vector_layout());

  REQUIRE (not fl1.is_tensor_layout());
  REQUIRE (not fl2.is_tensor_layout());
  REQUIRE (    fl3.is_tensor_layout());
  REQUIRE (not fl4.is_tensor_layout());
  REQUIRE (not fl5.is_tensor_layout());
  REQUIRE (    fl6.is_tensor_layout());

  REQUIRE (fl1.transpose()==fl1);
  REQUIRE (fl5.transpose()==fl7);

  REQUIRE (fl2.get_vector_tag()==CMP);
  REQUIRE (fl5.get_vector_tag()==CMP);
  REQUIRE (fl2.get_vector_component_idx()==1);
  REQUIRE (fl5.get_vector_component_idx()==1);
  REQUIRE (fl2.get_vector_dim()==1);
  REQUIRE (fl5.get_vector_dim()==2);

  REQUIRE (fl3.get_tensor_tags()==TVec{CMP,CMP});
  REQUIRE (fl6.get_tensor_components_ids()==IVec{1,2});
  REQUIRE (fl3.get_tensor_dims()==IVec{3,4});
  REQUIRE (fl6.get_tensor_dims()==IVec{5,6});
}

TEST_CASE("field_identifier", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags1 = {EL, LEV, CMP};
  std::vector<FieldTag> tags2 = {EL, CMP, LEV};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 4, 3};

  FieldIdentifier fid1 ("field_1", {tags1,dims1}, kg, "some_grid");
  FieldIdentifier fid2 ("field_1", {tags1,dims1}, kg, "some_grid");
  FieldIdentifier fid3 ("field_1", {tags1,dims2}, kg, "some_grid");
  FieldIdentifier fid4 ("field_2", {tags1,dims2}, kg, "some_grid");
  FieldIdentifier fid5 ("field_2", {tags2,dims2}, kg, "some_grid");
  FieldIdentifier fid6 ("field_2", {tags2,dims2}, m, "some_grid");
  FieldIdentifier fid7 ("field_2", {tags2,dims2}, m, "some_other_grid");

  REQUIRE (fid1==fid2);
  REQUIRE (fid2!=fid3);
  REQUIRE (fid3!=fid4);
  REQUIRE (fid4!=fid5);
  REQUIRE (fid5!=fid6);
  REQUIRE (fid6!=fid7);

  // Check that has_tag option works
  REQUIRE(fid1.get_layout().has_tag(CMP));
  REQUIRE(!fid1.get_layout().has_tag(GP));
}

TEST_CASE("field_tracking", "") {
  using namespace scream;

  FieldTracking track;
  util::TimeStamp time1(2021,10,12,17,8,10);
  util::TimeStamp time2(2021,10,12,17,8,20);
  track.update_time_stamp(time2);

  // Cannot rewind time (yet)
  REQUIRE_THROWS  (track.update_time_stamp(time1));
}

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

TEST_CASE("field_group") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FL  = FieldLayout;

  constexpr int ncols = 10;
  constexpr int ndims = 4;
  constexpr int nlevs = 8;

  FID fid ("V",FL({COL,CMP,LEV},{ncols,ndims,nlevs}),Units::nondimensional(),"the_grid");
  Field f (fid);
  f.allocate_view();

  FieldGroupInfo info("G");
  info.m_monolithic_allocation = true;
  std::vector<Field> f_i;
  for (int i=0; i<ndims; ++i) {
    f_i.push_back(f.get_component(i));
    info.m_fields_names.push_back(f_i[i].name());
    info.m_subview_dim = 1;
    info.m_subview_idx[f_i[i].name()] = i;
  }

  // Create group and set subfields
  FieldGroup g(info);
  g.m_monolithic_field = std::make_shared<Field>(f);
  for (int i=0; i<ndims; ++i) {
    g.m_individual_fields["G_"+std::to_string(i)] = std::make_shared<Field>(f_i[i]);
  }

  // Check const cloning
  auto cg= g.get_const();
  REQUIRE (cg.m_monolithic_field->is_read_only());
  REQUIRE (cg.m_individual_fields.size()==g.m_individual_fields.size());
  REQUIRE (*cg.m_info==*g.m_info);
  REQUIRE (cg.m_monolithic_field->get_internal_view_data<const Real>()==
            g.m_monolithic_field->get_internal_view_data<const Real>());
  for (int i=0; i<ndims; ++i) {
    const auto&  f =  *g.m_individual_fields.at("G_"+std::to_string(i));
    const auto& cf = *cg.m_individual_fields.at("G_"+std::to_string(i));
    REQUIRE ( f.get_internal_view_data<const Real>()==
             cf.get_internal_view_data<const Real>());
  }
}

TEST_CASE ("update") {
  using namespace scream;
  using namespace ekat::units;

  using namespace ShortFieldTagsNames;
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  auto engine = setup_random_test ();
  RPDF rpdf(0,1);
  IPDF ipdf(0,100);

  const int ncol = 2;
  const int ncmp = 3;
  const int nlev = 4;

  // Create field (if available, use packs, to ensure we don't print garbage)
  std::vector<FieldTag> tags = {COL, CMP, LEV};
  std::vector<int>      dims = {ncol,ncmp,nlev};

  FieldIdentifier fid_r ("fr", {tags,dims}, kg, "some_grid", DataType::RealType);
  FieldIdentifier fid_i ("fi", {tags,dims}, kg, "some_grid", DataType::IntType);
  Field f_real (fid_r);
  Field f_int  (fid_i);
  f_real.allocate_view();
  f_int.allocate_view();
  randomize (f_real,engine,rpdf);
  randomize (f_int, engine,ipdf);

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
    auto mask = f_int.clone("mask");
    compute_mask<Comparison::GT>(f1,0,mask);

    // Set f3=1 where mask=1
    f3.deep_copy(1,mask);

    auto one = f1.subfield(0,0).clone("one");
    auto zero = f1.subfield(0,0).clone("zero");
    one.deep_copy(1);
    zero.deep_copy(0);

    // Check
    for (int icol=0; icol<ncol; ++ icol) {
      auto f3i = f3.subfield(0,icol);
      if (icol % 2 == 0) {
        REQUIRE (views_are_equal(f3i,one));
      } else {
        REQUIRE (views_are_equal(f3i,zero));
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

TEST_CASE ("sync_subfields") {
  // This test is for previously incorrect behavior, where syncing a subfield
  // to host/device would deep copy the entire data view (including all entries of
  // the parent view). Here, if memory space is not shared between host and device,
  // syncing a subfield to host/device will not sync the data of the other subfields.

  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FL  = FieldLayout;

  constexpr int ncols = 10;
  constexpr int ndims = 4;
  constexpr int nlevs = 8;

  // Create field with (col, cmp, lev)
  FID fid ("V",FL({COL,CMP,LEV},{ncols,ndims,nlevs}),Units::nondimensional(),"the_grid",DataType::IntType);
  Field f (fid);
  f.allocate_view();

  // Store whether mem space for host and device are the same for testing subfield values
  const bool shared_mem_space = f.host_and_device_share_memory_space();

  // Deep copy all values to ndims on device and host
  f.deep_copy(ndims);
  f.sync_to_host();

  // Set subfield values to their index on device
  for (int c=0; c<ndims; ++c) {
    f.get_component(c).deep_copy(c);
  }

  // Sync only component 0 to host
  f.get_component(0).sync_to_host();

  // For components 1,...,ndims-1, if device and host do not share a
  // memory space, host values should be equal to ndims, else host
  // values should be equal to component index
  for (int c=1; c<ndims; ++c) {
    auto host_subview = f.get_component(c).get_view<int**, Host>();
    for (int idx=0; idx<ncols*nlevs; ++idx) {
      const int icol = idx/nlevs; const int ilev = idx%nlevs;
      if (shared_mem_space) REQUIRE(host_subview(icol, ilev) == c);
      else                  REQUIRE(host_subview(icol, ilev) == ndims);
    }
  }

  // Deep copy all values to ndims on device and host
  f.deep_copy(ndims);
  f.sync_to_host();

  // Set subfield values to their index on host
  for (int c=0; c<ndims; ++c) {
    f.get_component(c).deep_copy<Host>(c);
  }

  // Sync only component 0 to device
  f.get_component(0).sync_to_dev();

  // For components 1,...,ndims-1, if device and host do not share a
  // memory space, device values should be equal to ndims, else device
  // values should be equal to component index
  for (int c=1; c<ndims; ++c) {
    auto device_subview = f.get_component(c).get_view<int**, Device>();
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {ncols,nlevs}),
                         KOKKOS_LAMBDA (const int icol, const int ilev) {
      if (shared_mem_space) EKAT_KERNEL_ASSERT(device_subview(icol, ilev) == c);
      else                  EKAT_KERNEL_ASSERT(device_subview(icol, ilev) == ndims);
    });
  }
}

} // anonymous namespace
