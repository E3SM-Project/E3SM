#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_positivity_check.hpp"
#include "share/field/field_within_interval_check.hpp"
#include "share/field/field_monotonicity_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/user_provided_grids_manager.hpp"

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

#ifdef SCREAM_FORCE_RUN_FAIL
  REQUIRE(false); // force this test to fail
#endif

#ifdef SCREAM_FORCE_RUN_FPE_FAIL
  float foo = 42.0;
  float bar = foo / 0.0;
  std::cout << bar << std::endl;
#endif

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

  FieldTracking track("track");
  util::TimeStamp time1(0,0,0,10.0);
  util::TimeStamp time2(0,0,0,20.0);
  REQUIRE_NOTHROW (track.update_time_stamp(time2));

  REQUIRE_THROWS  (track.update_time_stamp(time1));
}

TEST_CASE("field", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;

  using P4 = ekat::Pack<Real,4>;
  using P8 = ekat::Pack<Real,8>;
  using P16 = ekat::Pack<Real,16>;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");

  // Check if we can extract a reshaped view
  SECTION ("reshape") {
    Field<Real> f1 (fid);

    // Should not be able to reshape before allocating
    REQUIRE_THROWS(f1.get_reshaped_view<Real*>());

    f1.allocate_view();

    // Reshape should work with both dynamic and static dims
    auto v1 = f1.get_reshaped_view<Real[3][24]>();
    auto v2 = f1.get_reshaped_view<Real**>();

    REQUIRE(v1.size()==v2.size());

    // But if wrong static length is used, we should throw
    REQUIRE_THROWS(f1.get_reshaped_view<Real[3][16]>());

    // Should not be able to reshape to this data type...
    REQUIRE_THROWS(f1.get_reshaped_view<P16**>());
    // But this should work
    REQUIRE_NOTHROW(f1.get_reshaped_view<P8**>());

    // Using packs (of allowable size) of different pack sizes
    // should lead to views with different extents.
    // Since there's no padding, their extent on last dimension
    // should be the phys dim divided by pack size.
    auto v3 = f1.get_reshaped_view<P8**>();
    auto v4 = f1.get_reshaped_view<P4**>();
    REQUIRE (v4.size() == 2*v3.size());
    REQUIRE (v4.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v3.extent_int(0) == fid.get_layout().dim(0));
    REQUIRE (v4.extent_int(1) == fid.get_layout().dim(1) / P4::n);
    REQUIRE (v3.extent_int(1) == fid.get_layout().dim(1) / P8::n);

    // The memory spans should be identical
    REQUIRE (v3.impl_map().memory_span()==v4.impl_map().memory_span());

    // Trying to reshape into something that the allocation cannot accommodate should throw
    REQUIRE_THROWS (f1.get_reshaped_view<P16***>());
  }

  SECTION ("compare") {

    Field<Real> f1(fid), f2(fid);
    f2.get_header().get_alloc_properties().request_allocation<P16>();
    f1.allocate_view();
    f2.allocate_view();

    auto v1 = f1.get_reshaped_view<Real**>();
    auto v2 = f2.get_reshaped_view<P8**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = i*dim1+j;

      int jpack = j / P8::n;
      int jvec = j % P8::n;
      v2(i,jpack)[jvec] = i*dim1+j;
    });
    Kokkos::fence();

    // The views were filled the same way, so they should test equal
    // NOTE: this cmp function only test the "actual" field, discarding padding.
    REQUIRE(views_are_equal(f1,f2));

    // Check self equivalence
    // get_const returns a copy of self, so equivalent (if already allocated)
    REQUIRE (f1.equivalent(f1.get_const()));
    REQUIRE (f1.equivalent(f1));
    // f1 and f2 have independent views, so they are not equivalent.
    REQUIRE (!f1.equivalent(f2));
  }

  // Check copy constructor
  SECTION ("copy ctor") {
    Field<Real> f1 (fid);

    f1.allocate_view();
    Kokkos::deep_copy(f1.get_view(),3.0);

    Field<const Real> f2 = f1;
    REQUIRE(f2.get_header_ptr()==f1.get_header_ptr());
    REQUIRE(f2.get_view()==f1.get_view());
    REQUIRE(f2.is_allocated());
    REQUIRE(views_are_equal(f1,f2));
  }

  SECTION ("deep_copy") {
    std::vector<FieldTag> t1 = {COL,CMP,LEV};
    std::vector<int> d1 = {3,2,24};

    FieldIdentifier fid1("vec_3d",{t1,d1},m/s,"some_grid");

    Field<Real> f1(fid1);
    f1.allocate_view();
    f1.deep_copy(1.0);
    f1.sync_to_host();
    auto v = f1.get_view<Host>();
    for (int i=0; i<fid1.get_layout().size(); ++i) {
      REQUIRE (v(i)==1.0);
    }
  }

  // Subfields
  SECTION ("subfield") {
    std::vector<FieldTag> t1 = {COL,CMP,CMP,LEV};
    std::vector<int> d1 = {3,10,2,24};

    FieldIdentifier fid1("4d",{t1,d1},m/s,"some_grid");

    Field<Real> f1(fid1);
    f1.allocate_view();
    auto v = f1.get_view();
    Kokkos::parallel_for(kt::RangePolicy(0,v.size()),
                         KOKKOS_LAMBDA(const int i) {
      v(i) = i;
    });
    auto vh = Kokkos::create_mirror_view(v);
    Kokkos::deep_copy(vh,v);

    const int idim = 1;
    const int ivar = 2;

    auto f2 = f1.subfield(idim,ivar);
    auto v4d = f1.get_reshaped_view<Real****>();
    auto v3d = f2.get_reshaped_view<Real***>();

    // Wrong rank for the subfield f2
    REQUIRE_THROWS(f2.get_reshaped_view<Real****>());

    auto v4d_h = f1.get_reshaped_view<Real****,Host>();
    auto v3d_h = f2.get_reshaped_view<Real***,Host>();
    for (int i=0; i<d1[0]; ++i)
      for (int j=0; j<d1[2]; ++j)
        for (int k=0; k<d1[3]; ++k) {
          REQUIRE (v4d_h(i,ivar,j,k)==v3d_h(i,j,k));
        }
  }

  SECTION ("vector_component") {
    std::vector<FieldTag> tags_2 = {COL,CMP,LEV};
    std::vector<int> dims_2 = {3,2,24};

    FieldIdentifier fid_2("vec_3d",{tags_2,dims_2},m/s,"some_grid");

    Field<Real> f_vec(fid_2);
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

    auto v = f_vec.get_reshaped_view<Real***,Host>();
    for (int col=0; col<3; ++col) {
      for (int lev=0; lev<24; ++lev) {
        REQUIRE (v(col,0,lev)==1.0);
        REQUIRE (v(col,1,lev)==2.0);
      }
    }
  }


  SECTION ("host_view") {
    Field<Real> f(fid);

    // Views not yet allocated
    REQUIRE_THROWS(f.get_view());
    REQUIRE_THROWS(f.get_view<Host>());
    REQUIRE_THROWS(f.sync_to_host());
    REQUIRE_THROWS(f.sync_to_dev());

    f.allocate_view();

    auto v = f.get_view();
    Kokkos::parallel_for(kt::RangePolicy(0,v.size()),
                         KOKKOS_LAMBDA(int i) {
      v(i) = i;
    });
    f.sync_to_host();

    // Get reshaped view on device, and manually create Host mirror
    auto v2d = f.get_reshaped_view<Real**>();
    auto v2d_hm = Kokkos::create_mirror_view(v2d);
    Kokkos::deep_copy(v2d_hm,v2d);

    // Get reshaped view straight on Host
    auto v2dh = f.get_reshaped_view<Real**,Host>();

    // The two should match
    for (int i=0; i<dims[0]; ++i) {
      for (int j=0; j<dims[1]; ++j) {
        REQUIRE (v2dh(i,j) == v2d_hm(i,j) );
      }
    }
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

  std::vector<FieldTag> tags1 = {COL,LEV};
  std::vector<FieldTag> tags2 = {EL,GP,GP};

  std::vector<int> dims1 = {ncols,nlevs};
  std::vector<int> dims2 = {2,4,4};
  std::vector<int> dims3 = {ncols,nlevs+1};

  const auto km = 1000*m;

  FID fid1("field_1", {tags1, dims1},  m/s, "phys");
  FID fid2("field_2", {tags1, dims1},  m/s, "phys");
  FID fid3("field_3", {tags1, dims1},  m/s, "phys");

  FID bad1("field_1", {tags2, dims2},  m/s, "dyn");  // Bad grid
  FID bad2("field_2", {tags1, dims1}, km/s, "phys"); // Bad units
  FID bad3("field_2", {tags1, dims3},  m/s, "phys"); // Bad layout

  ekat::Comm comm(MPI_COMM_WORLD);
  auto pg = create_point_grid("phys",ncols*comm.size(),nlevs,comm);
  FieldManager<Real> field_mgr(pg);

  // Should not be able to register fields yet
  REQUIRE_THROWS(field_mgr.register_field(FR(fid1)));

  field_mgr.registration_begins();

  // === Valid registration calls === //
  field_mgr.register_field(FR(fid1,"group_1",Pack::n));
  field_mgr.register_field(FR{fid2,"group_2",16});
  field_mgr.register_field(FR{fid3,"group_4"});
  field_mgr.register_field(FR{fid3,SL{"group_1","group_2","group_3"}});
  field_mgr.register_field(FR{fid2,"group_4"});

  // === Invalid registration calls === //
  REQUIRE_THROWS(field_mgr.register_field(FR{bad1}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));
  field_mgr.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(field_mgr.register_field(FR{fid1}));

  // Check registration is indeed closed
  REQUIRE (field_mgr.repository_state()==RepoState::Closed);
  REQUIRE (field_mgr.size()==3);

  // Get all fields
  auto f1 = field_mgr.get_field(fid1.name());
  auto f2 = field_mgr.get_field(fid2.name());
  auto f3 = field_mgr.get_field(fid3.name());
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

  FieldManager<Real> field_mgr(pg);
  field_mgr.registration_begins();
  field_mgr.register_field(FR{qv_id,"tracers"});
  field_mgr.register_field(FR{qc_id,"tracers"});
  field_mgr.register_field(FR{qr_id,"tracers"});
  field_mgr.register_group(GroupRequest("tracers",grid_name,Bundling::Required));
  field_mgr.registration_ends();

  auto qv = field_mgr.get_field(qv_id.name());
  auto qc = field_mgr.get_field(qc_id.name());
  auto qr = field_mgr.get_field(qr_id.name());

  auto group = field_mgr.get_field_group("tracers");
  // The field_mgr should have allocated the group bundled
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
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);

  ekat::genRandArray(Q.get_view(),engine,pdf);

  // Check that the same values are in all q's
  Q.sync_to_host();
  auto Qh = Q.get_reshaped_view<Real***,Host>();
  auto qvh = qv.get_reshaped_view<Real**,Host>();
  auto qch = qc.get_reshaped_view<Real**,Host>();
  auto qrh = qr.get_reshaped_view<Real**,Host>();

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

  FieldRequest a_req(a_id,"group1");
  FieldRequest b_req(b_id,SL{"group1"});
  FieldRequest c_req(c_id,SL{"group1","group2"});
  FieldRequest d_req(d_id,SL{"group1","group3"});
  FieldRequest e_req(e_id,SL{"group1","group2"});
  FieldRequest f_req(f_id,SL{"group1","group3"});

  GroupRequest g1_req ("group1",grid_name,Bundling::Required);
  GroupRequest g2_req ("group2",grid_name,Bundling::Required);
  // group3 = group3 + group2
  GroupRequest g3_req ("group3",grid_name,4,Bundling::Required,&g2_req,Relationship::Child);
  // group4 = group2
  GroupRequest g4_req ("group4",grid_name,4,Bundling::Required,&g2_req,Relationship::Alias);
  // group5 = group1 - {c,d}
  GroupRequest g5_req ("group5",grid_name,4,Bundling::Preferred,&g1_req,Relationship::Parent,SL{"c","d"});

  // The above group specs give the following groups:
  // g1: [a,b,c,d,e,f]
  // g2: [c,e]
  // g3: [d,f,c,e]
  // g4: [c,e]
  // g5: [a,b,e,f]
  // The bundling requests can be accommodated for g1-g4, but not g5.
  // But g5 request is only 'Preferred', so the FM won't error out.
  // The order of fields in the 'encompassing' group is {[c,e][d,f][a,b]},
  // where [c,e] means that the order of those two fields can be anything.
  // The 'block'-reverse of that list is also possible: {[a,b][d,f][c,e]}

  FieldManager<Real> field_mgr(pg);
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
  if (f1=="a" || f1=="b") {
    REQUIRE ( ((f1=="a" && f2=="b") || (f1=="b" && f2=="a")) );
    REQUIRE ( ((f3=="d" && f4=="f") || (f3=="f" && f4=="d")) );
    REQUIRE ( ((f5=="c" && f6=="e") || (f5=="e" && f6=="c")) );
  } else {
    REQUIRE ( ((f1=="c" && f2=="e") || (f1=="e" && f2=="c")) );
    REQUIRE ( ((f3=="d" && f4=="f") || (f3=="f" && f4=="d")) );
    REQUIRE ( ((f5=="a" && f6=="b") || (f5=="b" && f6=="a")) );
  }
}

TEST_CASE("field_property_check", "") {

  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags = {EL, GP, LEV};
  std::vector<int> dims = {2, 3, 12};

  FieldIdentifier fid ("field_1",{tags,dims}, m/s,"some_grid");

  // Check positivity.
  SECTION ("field_positivity_check") {
    Field<Real> f1(fid);
    auto positivity_check = std::make_shared<FieldPositivityCheck<Real> >();
    REQUIRE(not positivity_check->can_repair());
    f1.add_property_check(positivity_check);
    f1.allocate_view();

    // Assign positive values to the field and make sure it passes our test for
    // positivity.
    auto f1_view = f1.get_view();
    auto host_view = Kokkos::create_mirror_view(f1_view);
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = i+1;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(iter->check(f1));
    }

    // Assign non-positive values to the field and make sure it fails the check.
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = -i;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(not iter->check(f1));
    }
  }

  // Check positivity with repairs.
  SECTION ("field_positivity_check_with_repairs") {
    Field<Real> f1(fid);
    auto positivity_check = std::make_shared<FieldPositivityCheck<Real> >(1);
    REQUIRE(positivity_check->can_repair());
    f1.add_property_check(positivity_check);
    f1.allocate_view();

    // Assign non-positive values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    auto f1_view = f1.get_view();
    auto host_view = Kokkos::create_mirror_view(f1_view);
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = -i;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(not iter->check(f1));
      iter->repair(f1);
      REQUIRE(iter->check(f1));
    }
  }

  // Check that the values of a field lie within an interval.
  SECTION ("field_within_interval_check") {
    Field<Real> f1(fid);
    auto interval_check = std::make_shared<FieldWithinIntervalCheck<Real> >(0, 100);
    REQUIRE(interval_check->can_repair());
    f1.add_property_check(interval_check);
    f1.allocate_view();

    // Assign positive values to the field and make sure it passes our test for
    // positivity.
    auto f1_view = f1.get_view();
    auto host_view = Kokkos::create_mirror_view(f1_view);
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = i;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(iter->check(f1));
    }

    // Assign non-positive values to the field, make sure it fails the check,
    // and then repair the field so it passes.
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = -i;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(not iter->check(f1));
      iter->repair(f1);
      REQUIRE(iter->check(f1));
    }
  }

  // Check monotonicity.
  SECTION ("field_monotonicity_check") {
    Field<Real> f1(fid);
    auto mono_check = std::make_shared<FieldMonotonicityCheck<Real> >();
    REQUIRE(not mono_check->can_repair());
    f1.add_property_check(mono_check);
    f1.allocate_view();

    // Assign monotonically-increasing values to the field and make sure it
    // passes our test for positivity.
    auto f1_view = f1.get_view();
    auto host_view = Kokkos::create_mirror_view(f1_view);
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = i+1;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(iter->check(f1));
    }

    // Assign monotonically-decreasing values to the field and make sure it
    // also passes the check.
    for (int i = 0; i < host_view.extent_int(0); ++i) {
      host_view(i) = -i;
    }
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(iter->check(f1));
    }

    // Write a positive value to the middle of the array that causes the
    // monotonicity check to fail.
    host_view(host_view.extent(0)/2) = 1;
    Kokkos::deep_copy(f1_view, host_view);
    for (auto iter = f1.property_check_begin(); iter != f1.property_check_end(); iter++) {
      REQUIRE(not iter->check(f1));
      REQUIRE_THROWS(iter->repair(f1)); // we can't repair it, either
    }
  }
}

} // anonymous namespace
