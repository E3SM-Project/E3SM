#include <catch2/catch.hpp>
#include <numeric>

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "share/grid/point_grid.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("field_layout") {
  using namespace scream;
  using namespace ShortFieldTagsNames;

  FieldLayout l({EL,GP,GP});

  // Should not be able to set a dimensions vector of wrong rank
  REQUIRE_THROWS(l.set_dimensions({1,2}));

  l.set_dimensions({1,2,3});

  // Should not be able to reset the dimensions once they are set
  REQUIRE_THROWS(l.set_dimensions({1,2,3}));
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
  REQUIRE_NOTHROW (track.update_time_stamp(time2));

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
    REQUIRE_NOTHROW(f1.get_view<P8**>());

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
  }

  SECTION ("equivalent") {
    Field f1 (fid), f2(fid);
    f1.allocate_view();
    f2.allocate_view();

    // Check self equivalence
    // get_const returns a copy of self, so equivalent (if already allocated)
    REQUIRE (f1.equivalent(f1.get_const()));
    REQUIRE (f1.equivalent(f1));
    // f1 and f2 have independent views, so they are not equivalent.
    REQUIRE (!f1.equivalent(f2));
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
    f2.deep_copy<Real>(0.0);
    REQUIRE (field_max<Real>(f2)==0.0);
    REQUIRE (field_min<Real>(f2)==0.0);
    REQUIRE (field_max<Real>(f1)==3.0);
    REQUIRE (field_min<Real>(f1)==3.0);
  }

  SECTION ("deep_copy") {
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

  // Subfields
  SECTION ("subfield") {
    std::vector<FieldTag> t1 = {COL,CMP,CMP,LEV};
    std::vector<int> d1 = {3,10,2,24};

    FieldIdentifier fid1("4d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1,engine,pdf);

    const int idim = 1;
    const int ivar = 2;

    auto f2 = f1.subfield(idim,ivar);

    // Wrong rank for the subfield f2
    REQUIRE_THROWS(f2.get_view<Real****>());

    auto v4d_h = f1.get_view<Real****,Host>();
    auto v3d_h = f2.get_view<Real***,Host>();
    for (int i=0; i<d1[0]; ++i)
      for (int j=0; j<d1[2]; ++j)
        for (int k=0; k<d1[3]; ++k) {
          REQUIRE (v4d_h(i,ivar,j,k)==v3d_h(i,j,k));
        }
  }

  // Dynamic Subfields
  SECTION ("dynamic_subfield") {
    const int vec_dim = 10;
    std::vector<FieldTag> t1 = {COL,CMP,CMP,LEV};
    std::vector<int> d1 = {3,vec_dim,2,24};

    FieldIdentifier fid1("4d",{t1,d1},m/s,"some_grid");

    Field f1(fid1);
    f1.allocate_view();
    randomize(f1,engine,pdf);

    const int idim = 1;
    const int ivar = 0;

    auto f2 = f1.subfield(idim,ivar,/* dynamic = */ true);

    // Cannot reset subview idx of non-subfield fields
    REQUIRE_THROWS(f1.get_header().get_alloc_properties().reset_subview_idx(0));

    // subview idx out of bounds
    auto& f2_ap = f2.get_header().get_alloc_properties();
    REQUIRE_THROWS(f2_ap.reset_subview_idx(-1));
    REQUIRE_THROWS(f2_ap.reset_subview_idx(vec_dim));

    // Fill f1 with random numbers, and verify corresponding subviews get same values
    randomize(f1,engine,pdf);

    for (int ivar_dyn=0; ivar_dyn<vec_dim; ++ivar_dyn) {
      // Reset slice idx
      f2_ap.reset_subview_idx(ivar_dyn);
      REQUIRE(f2_ap.get_subview_info().slice_idx==ivar_dyn);

      auto v4d_h = f1.get_view<Real****,Host>();
      auto v3d_h = f2.get_view<Real***,Host>();
      for (int i=0; i<d1[0]; ++i)
        for (int j=0; j<d1[2]; ++j)
          for (int k=0; k<d1[3]; ++k) {
            REQUIRE (v4d_h(i,ivar_dyn,j,k)==v3d_h(i,j,k));
          }
    }
  }

  SECTION ("vector_component") {
    std::vector<FieldTag> tags_2 = {COL,CMP,LEV};
    std::vector<int> dims_2 = {3,2,24};

    FieldIdentifier fid_2("vec_3d",{tags_2,dims_2},m/s,"some_grid");

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

    auto v = f_vec.get_view<Real***,Host>();
    for (int col=0; col<3; ++col) {
      for (int lev=0; lev<24; ++lev) {
        REQUIRE (v(col,0,lev)==1.0);
        REQUIRE (v(col,1,lev)==2.0);
      }
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
  info.m_bundled = true;
  std::vector<Field> f_i;
  for (int i=0; i<ndims; ++i) {
    f_i.push_back(f.get_component(i));
    info.m_fields_names.push_back(f_i[i].name());
    info.m_subview_dim = 1;
    info.m_subview_idx[f_i[i].name()] = i;
  }

  // Create group and set subfields
  FieldGroup g(info);
  g.m_bundle = std::make_shared<Field>(f);
  for (int i=0; i<ndims; ++i) {
    g.m_fields["G_"+std::to_string(i)] = std::make_shared<Field>(f_i[i]);
  }

  // Check const cloning
  auto cg= g.get_const();
  REQUIRE (cg.m_bundle->is_read_only());
  REQUIRE (cg.m_fields.size()==g.m_fields.size());
  REQUIRE (*cg.m_info==*g.m_info);
  REQUIRE (cg.m_bundle->get_internal_view_data<const Real>()==
            g.m_bundle->get_internal_view_data<const Real>());
  for (int i=0; i<ndims; ++i) {
    const auto&  f =  *g.m_fields.at("G_"+std::to_string(i));
    const auto& cf = *cg.m_fields.at("G_"+std::to_string(i));
    REQUIRE ( f.get_internal_view_data<const Real>()==
             cf.get_internal_view_data<const Real>());
  }
}

TEST_CASE("field_mgr", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FR  = FieldRequest;
  using SL  = std::list<std::string>;
  using Pack = ekat::Pack<Real,8>;

  const int ncols = 4;
  const int nlevs = 7;
  const int subview_dim = 1;
  const int subview_slice = 0;

  std::vector<FieldTag> tags1 = {COL,LEV};
  std::vector<FieldTag> tags2 = {EL,GP,GP};
  std::vector<FieldTag> tags3 = {COL,CMP,LEV};

  std::vector<int> dims1 = {ncols,nlevs};
  std::vector<int> dims2 = {2,4,4};
  std::vector<int> dims3 = {ncols,nlevs+1};
  std::vector<int> dims4 = {ncols,10,nlevs};

  const auto km = 1000*m;

  FID fid1("field_1", {tags1, dims1},  m/s, "phys");
  FID fid2("field_2", {tags1, dims1},  m/s, "phys");
  FID fid3("field_3", {tags1, dims1},  m/s, "phys");
  FID fid4("field_4", {tags3, dims4},  m/s, "phys");

  FID bad1("field_1", {tags2, dims2},  m/s, "dyn");  // Bad grid
  FID bad2("field_2", {tags1, dims1}, km/s, "phys"); // Bad units
  FID bad3("field_2", {tags1, dims3},  m/s, "phys"); // Bad layout

  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid("phys",ncols*comm.size(),nlevs,comm);
  FieldManager field_mgr(pg);

  // Should not be able to register fields yet
  REQUIRE_THROWS(field_mgr.register_field(FR(fid1)));

  field_mgr.registration_begins();

  // === Valid registration calls === //
  field_mgr.register_field(FR(fid1,"group_1",Pack::n));
  field_mgr.register_field(FR{fid2,"group_2",16});
  field_mgr.register_field(FR{fid3,"group_4"});
  field_mgr.register_field(FR{fid3,SL{"group_1","group_2","group_3"}});
  field_mgr.register_field(FR{fid2,"group_4"});
  field_mgr.register_field(FR{fid4});

  // === Invalid registration calls === //
  REQUIRE_THROWS(field_mgr.register_field(FR{bad1}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));

  // Cannot add external fields while registration is happening
  REQUIRE_THROWS(field_mgr.add_field(Field()));

  field_mgr.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(field_mgr.register_field(FR{fid1}));

  // Check registration is indeed closed
  REQUIRE (field_mgr.repository_state()==RepoState::Closed);
  REQUIRE (field_mgr.size()==4);

  // Get all fields
  auto f1 = field_mgr.get_field(fid1.name());
  auto f2 = field_mgr.get_field(fid2.name());
  auto f3 = field_mgr.get_field(fid3.name());
  auto f4 = field_mgr.get_field(fid4.name());
  REQUIRE_THROWS(field_mgr.get_field("bad")); // Not in the field_mgr
  REQUIRE(f1.get_header().get_identifier()==fid1);

  // Check that the groups names are in the header. While at it, make sure that case insensitive works fine.
  auto has_group = [](const ekat::WeakPtrSet<const FieldGroupInfo>& groups,
                      const std::string& name)->bool {
    for (auto it : groups) {
      if (it.lock()->m_group_name==name) {
        return true;
      }
    }
    return false;
  };
  REQUIRE (has_group(f1.get_header().get_tracking().get_groups_info(),"gRouP_1"));
  REQUIRE (has_group(f2.get_header().get_tracking().get_groups_info(),"Group_2"));
  REQUIRE (has_group(f2.get_header().get_tracking().get_groups_info(),"Group_4"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_1"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_2"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_3"));
  REQUIRE (has_group(f3.get_header().get_tracking().get_groups_info(),"Group_4"));

  // Check that the groups in the field_mgr contain the correct fields
  REQUIRE (field_mgr.get_groups_info().count("GROUP_1")==1);
  REQUIRE (field_mgr.get_groups_info().count("GRoup_2")==1);
  REQUIRE (field_mgr.get_groups_info().count("group_3")==1);
  REQUIRE (field_mgr.get_groups_info().count("groUP_4")==1);
  REQUIRE (field_mgr.get_groups_info().count("group_5")==0);
  REQUIRE (field_mgr.get_groups_info().at("group_2")->m_fields_names.size()==2);

  auto g1 = field_mgr.get_groups_info().at("group_1");
  auto g2 = field_mgr.get_groups_info().at("group_2");
  auto g3 = field_mgr.get_groups_info().at("group_3");
  auto g4 = field_mgr.get_groups_info().at("group_4");
  REQUIRE (ekat::contains(g1->m_fields_names,"Field_1"));
  REQUIRE (ekat::contains(g1->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g2->m_fields_names,"Field_2"));
  REQUIRE (ekat::contains(g2->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g3->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(g4->m_fields_names,"Field_2"));
  REQUIRE (ekat::contains(g4->m_fields_names,"Field_3"));

  // Check alloc props for f1 and f2 (which requested pack size > 1)
  auto f1_padding = f1.get_header().get_alloc_properties().get_padding();
  auto f2_padding = f2.get_header().get_alloc_properties().get_padding();

  REQUIRE (f1_padding==ekat::PackInfo<Pack::n>::padding(nlevs));
  REQUIRE (f2_padding==ekat::PackInfo<16>::padding(nlevs));

  // Try to subview a field and set the subfield back in the FM
  field_mgr.add_field(f4.subfield("field_4_sf",subview_dim,subview_slice,true));
  auto f4_sf = field_mgr.get_field("field_4_sf");
  REQUIRE (field_mgr.size()==5);
  REQUIRE_THROWS (field_mgr.add_field(Field())); // Not allocated
  REQUIRE_THROWS (field_mgr.add_field(f4_sf)); // Cannot have duplicates

  // Verify f5 is a subfield of f4
  auto f4_sf_ap = f4_sf.get_header().get_alloc_properties();
  REQUIRE (f4_sf_ap.is_subfield());
  REQUIRE (f4_sf_ap.is_dynamic_subfield());
  REQUIRE (f4_sf_ap.get_subview_info().dim_idx==subview_dim);
  REQUIRE (f4_sf_ap.get_subview_info().slice_idx==subview_slice);

  // Fill f4_sf with random numbers, and verify corresponding subview of f4 gets same values.
  auto engine = setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);
  randomize(f4_sf,engine,pdf);
  REQUIRE (views_are_equal(f4_sf,f4.get_component(subview_slice)));
}

TEST_CASE("tracers_bundle", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR  = FieldRequest;

  const int ncols = 4;
  const int nlevs = 7;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {ncols,nlevs};

  const auto nondim = Units::nondimensional();

  const std::string grid_name = "physics";
  FieldIdentifier qv_id("qv", {tags, dims}, nondim, grid_name);
  FieldIdentifier qc_id("qc", {tags, dims}, nondim, grid_name);
  FieldIdentifier qr_id("qr", {tags, dims}, nondim, grid_name);

  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid(grid_name,ncols*comm.size(),nlevs,comm);

  FieldManager field_mgr(pg);
  field_mgr.registration_begins();
  field_mgr.register_field(FR{qv_id,"tracers"});
  field_mgr.register_field(FR{qc_id,"tracers"});
  field_mgr.register_field(FR{qr_id,"tracers"});
  field_mgr.register_group(GroupRequest("tracers",grid_name,Bundling::Required));
  field_mgr.registration_ends();

  auto qv = field_mgr.get_field(qv_id.name());
  auto qc = field_mgr.get_field(qc_id.name());
  auto qr = field_mgr.get_field(qr_id.name());

  // The field_mgr should have allocated the group bundled
  auto group = field_mgr.get_field_group("tracers");
  REQUIRE (group.m_info->m_bundled);

  const auto& Q_name = group.m_bundle->get_header().get_identifier().name();
  auto Q = field_mgr.get_field(Q_name);

  // The bundled field in the group should match the field we get from the field_mgr
  REQUIRE (Q.equivalent(*group.m_bundle));

  // Check that Q is set as parent for all q's.
  auto qvp = qv.get_header().get_parent().lock();
  auto qcp = qc.get_header().get_parent().lock();
  auto qrp = qr.get_header().get_parent().lock();
  REQUIRE ((qvp!=nullptr && qvp.get()==&Q.get_header()));
  REQUIRE ((qcp!=nullptr && qvp.get()==&Q.get_header()));
  REQUIRE ((qrp!=nullptr && qvp.get()==&Q.get_header()));

  // The indices used for each q to subview Q
  int idx_v, idx_c, idx_r;

  // The idx must be stored
  REQUIRE_NOTHROW (idx_v = group.m_info->m_subview_idx.at("qv"));
  REQUIRE_NOTHROW (idx_c = group.m_info->m_subview_idx.at("qc"));
  REQUIRE_NOTHROW (idx_r = group.m_info->m_subview_idx.at("qr"));

  // All idx must be in [0,2] and must be different
  REQUIRE ((idx_v>=0 && idx_v<3 &&
            idx_c>=0 && idx_c<3 &&
            idx_r>=0 && idx_r<3));
  REQUIRE ((idx_v!=idx_c && idx_v!=idx_r && idx_c!=idx_r));

  // Now fill Q with random values
  auto engine = setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);

  randomize(Q,engine,pdf);

  // Check that the same values are in all q's
  Q.sync_to_host();
  auto Qh = Q.get_view<Real***,Host>();
  auto qvh = qv.get_view<Real**,Host>();
  auto qch = qc.get_view<Real**,Host>();
  auto qrh = qr.get_view<Real**,Host>();

  for (int icol=0; icol<ncols; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      REQUIRE (Qh(icol,idx_v,ilev)==qvh(icol,ilev));
      REQUIRE (Qh(icol,idx_c,ilev)==qch(icol,ilev));
      REQUIRE (Qh(icol,idx_r,ilev)==qrh(icol,ilev));
    }
  }

  // Check that the field ptrs stored in the group are the same as the q
  auto qv_ptr = group.m_fields.at("qv");
  auto qc_ptr = group.m_fields.at("qc");
  auto qr_ptr = group.m_fields.at("qr");

  REQUIRE (qv_ptr->equivalent(qv));
  REQUIRE (qc_ptr->equivalent(qc));
  REQUIRE (qr_ptr->equivalent(qr));
}

TEST_CASE("multiple_bundles") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using SL = std::list<std::string>;

  const int ncols = 4;
  const int nlevs = 7;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {ncols,nlevs};

  const auto nondim = Units::nondimensional();

  const std::string grid_name = "physics";
  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid(grid_name,ncols*comm.size(),nlevs,comm);

  FieldIdentifier a_id("a", {tags, dims}, nondim, grid_name);
  FieldIdentifier b_id("b", {tags, dims}, nondim, grid_name);
  FieldIdentifier c_id("c", {tags, dims}, nondim, grid_name);
  FieldIdentifier d_id("d", {tags, dims}, nondim, grid_name);
  FieldIdentifier e_id("e", {tags, dims}, nondim, grid_name);
  FieldIdentifier f_id("f", {tags, dims}, nondim, grid_name);

  FieldRequest a_req(a_id,SL{"group1","group3"});
  FieldRequest b_req(b_id,SL{"group1"});
  FieldRequest c_req(c_id,SL{"group1","group2"});
  FieldRequest d_req(d_id,SL{"group1","group3"});
  FieldRequest e_req(e_id,SL{"group1","group2"});
  FieldRequest f_req(f_id,SL{"group1"});

  GroupRequest g1_req ("group1",grid_name,Bundling::Required);
  GroupRequest g2_req ("group2",grid_name,Bundling::Required);
  // Include all group2 in group3
  GroupRequest g3_req ("group3",grid_name,4,Bundling::Required,DerivationType::Superset,g2_req.name,g2_req.grid);
  // Create group4 as a copy of group2
  GroupRequest g4_req ("group4",grid_name,4,Bundling::Required,DerivationType::Copy,g2_req.name,g2_req.grid);
  // Extend group5 by adding all fields in group1 *except* 'c' and 'd'.
  GroupRequest g5_req ("group5",grid_name,4,Bundling::Preferred,DerivationType::Subset,g1_req.name,g1_req.grid,SL{"c","d"});

  // The above group specs should give the following groups:
  // g1: [a,b,c,d,e,f]
  // g2: [c,e]
  // g3: [a,c,d,e]
  // g4: [c,e]
  // g5: [a,b,e,f]
  // The bundling requests can be accommodated for g1,g2,g3,g4, but not g5.
  // But g5 request is only 'Preferred', so the FM won't error out.
  // The order of fields in the 'encompassing' group is {[c,e],[a,d],[b,f]},
  // where [f1,..,fn] means that the order of those two fields can be anything.
  // The 'block'-reverse of that list is also possible: {[b,f],[a,d],[c,e]}

  FieldManager field_mgr(pg);
  field_mgr.registration_begins();

  // Register single fields
  field_mgr.register_field(a_req);
  field_mgr.register_field(b_req);
  field_mgr.register_field(c_req);
  field_mgr.register_field(d_req);
  field_mgr.register_field(e_req);
  field_mgr.register_field(f_req);

  // Register groups
  field_mgr.register_group(g1_req);
  field_mgr.register_group(g2_req);
  field_mgr.register_group(g3_req);
  field_mgr.register_group(g4_req);
  field_mgr.register_group(g5_req);

  field_mgr.registration_ends();

  auto g1 = field_mgr.get_field_group(g1_req.name);
  auto g2 = field_mgr.get_field_group(g2_req.name);
  auto g3 = field_mgr.get_field_group(g3_req.name);
  auto g4 = field_mgr.get_field_group(g4_req.name);
  auto g5 = field_mgr.get_field_group(g5_req.name);

  // First 4 groups should be bundled
  REQUIRE (g1.m_info->m_bundled);
  REQUIRE (g2.m_info->m_bundled);
  REQUIRE (g3.m_info->m_bundled);
  REQUIRE (g4.m_info->m_bundled);
  REQUIRE (not g5.m_info->m_bundled);

  // Check that the order of fields in g1 is the expected one
  const auto& fnames = g1.m_info->m_fields_names;
  const auto& f1 = *std::next(fnames.begin(),0);
  const auto& f2 = *std::next(fnames.begin(),1);
  const auto& f3 = *std::next(fnames.begin(),2);
  const auto& f4 = *std::next(fnames.begin(),3);
  const auto& f5 = *std::next(fnames.begin(),4);
  const auto& f6 = *std::next(fnames.begin(),5);
  if (f1=="b" || f1=="f") {
    REQUIRE ( ((f1=="b" && f2=="f") || (f1=="f" && f2=="b")) );
    REQUIRE ( ((f3=="a" && f4=="d") || (f3=="d" && f4=="a")) );
    REQUIRE ( ((f5=="c" && f6=="e") || (f5=="e" && f6=="c")) );
  } else {
    REQUIRE ( ((f1=="c" && f2=="e") || (f1=="e" && f2=="c")) );
    REQUIRE ( ((f3=="a" && f4=="d") || (f3=="d" && f4=="a")) );
    REQUIRE ( ((f5=="b" && f6=="f") || (f5=="f" && f6=="b")) );
  }
}

} // anonymous namespace
