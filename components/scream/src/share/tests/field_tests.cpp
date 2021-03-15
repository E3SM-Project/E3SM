#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"
#include "share/field/field_positivity_check.hpp"
#include "share/field/field_within_interval_check.hpp"
#include "share/field/field_monotonicity_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/user_provided_grids_manager.hpp"

#include "ekat/ekat_pack.hpp"
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

  // Subfields
  SECTION ("subfield") {
    std::vector<FieldTag> t1 = {COL,VAR,CMP,LEV};
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

TEST_CASE("field_repo", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::vector<FieldTag> tags1 = {EL, GP, GP};
  std::vector<FieldTag> tags2 = {COL};
  std::vector<FieldTag> tags3 = {ILEV};
  std::vector<FieldTag> tags4 = {COL,LEV};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};
  std::vector<int> dims3 = {13};
  std::vector<int> dims4 = {6};
  std::vector<int> dims5 = {13,6};

  const auto km = 1000*m;

  FieldIdentifier fid1("field_1", {tags1, dims1},  m/s, "grid_1");
  FieldIdentifier fid2("field_2", {tags1, dims2},  m/s, "grid_1");
  FieldIdentifier fid3("field_2", {tags2, dims3},  m/s, "grid_2");
  FieldIdentifier fid4("field_2", {tags2, dims3}, km/s, "grid_1");
  FieldIdentifier fid5("field_3", {tags2, dims3}, km/s, "grid_1");
  FieldIdentifier fid6("field_4", {tags2, dims3}, km/s, "grid_1");
  FieldIdentifier fid7("field_5", {tags2, dims3}, km/s, "grid_3");
  FieldIdentifier fid8("field_5", {tags3, dims4}, km/s, "grid_3");
  FieldIdentifier fid9("field_packed", {tags4,dims5}, km/s, "grid_3");

  FieldRepository<Real> repo;

  // Should not be able to register fields yet
  REQUIRE_THROWS(repo.register_field(fid1,"group_1"));

  repo.registration_begins();
  repo.register_field(fid1,"group_1");
  repo.register_field(fid2,"group_2");
  repo.register_field(fid3,"group_2");
  // Test that you can assign more than one group to a field
  repo.register_field(fid5,{"group_3","group_5"});
  repo.register_field(fid6,{"group_4","group_5","group_6"});
  // Test for same field and grid name, different layout
  repo.register_field(fid7,"group_7");
  repo.register_field(fid8,"group_7");
  // Test for packed field
  using Pack         = ekat::Pack<Real,8>;
  repo.register_field<Pack>(fid9);
  // Should not be able to register the same field name with two different units
  REQUIRE_THROWS(repo.register_field(fid4));
  repo.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(repo.register_field(fid1,"group_1"));

  // Check registration is indeed closed
  REQUIRE (repo.repository_state()==RepoState::Closed);
  REQUIRE (repo.size()==6);
  REQUIRE (repo.internal_size()==8);

  auto f1 = repo.get_field(fid1);
  auto f2 = repo.get_field(fid2);
  auto f5 = repo.get_field(fid5);
  auto f6 = repo.get_field(fid6);

  // Check that get_field with a field name and grid name as arguments returns the appropriate field
  // Using grid_1 should return fid2 field, using grid_2 should return fid3 field, using grid_3 should throw an error.
  // Retrieving field 5 on grid 3 should return an error since it has been defined on grid 3 with two different layouts.
  auto f7 = repo.get_field("field_2","grid_1");
  auto f8 = repo.get_field("field_2","grid_2");
  REQUIRE_THROWS( repo.get_field("field_5","grid_3") );
  REQUIRE_THROWS( repo.get_field("field_2","grid_3") );
  REQUIRE(f7.get_header().get_identifier()==fid2);
  REQUIRE(f7.get_header().get_identifier()!=fid3);
  REQUIRE(f8.get_header().get_identifier()!=fid2);
  REQUIRE(f8.get_header().get_identifier()==fid3);

  // Check the two fields identifiers are indeed different
  REQUIRE (f1.get_header().get_identifier()!=f2.get_header().get_identifier());

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
  REQUIRE (has_group(f5.get_header().get_tracking().get_groups_info(),"Group_3"));
  REQUIRE (has_group(f5.get_header().get_tracking().get_groups_info(),"Group_5"));
  REQUIRE (has_group(f6.get_header().get_tracking().get_groups_info(),"Group_4"));
  REQUIRE (has_group(f6.get_header().get_tracking().get_groups_info(),"Group_5"));
  REQUIRE (has_group(f6.get_header().get_tracking().get_groups_info(),"Group_6"));

  // Check that the groups in the repo contain the correct fields
  REQUIRE (repo.get_groups_info().count("GROUP_1")==1);
  REQUIRE (repo.get_groups_info().count("GRoup_2")==1);
  REQUIRE (repo.get_groups_info().count("group_3")==1);
  REQUIRE (repo.get_groups_info().count("groUP_4")==1);
  REQUIRE (repo.get_groups_info().count("group_5")==1);
  REQUIRE (repo.get_groups_info().at("group_5")->m_fields_names.size()==2);
  REQUIRE (repo.get_groups_info().count("group_6")==1);
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_1")->m_fields_names,"Field_1"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_2")->m_fields_names,"Field_2"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_3")->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_4")->m_fields_names,"Field_4"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_5")->m_fields_names,"Field_3"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_5")->m_fields_names,"Field_4"));
  REQUIRE (ekat::contains(repo.get_groups_info().at("group_6")->m_fields_names,"Field_4"));

  // Check that get_padding returns the appropriate value
  auto f9 = repo.get_field("field_packed","grid_3");
  REQUIRE (f9.get_header().get_alloc_properties().get_padding()==2);
}

TEST_CASE("tracers_bundle", "") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

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
  auto gm = create_user_provided_grids_manager(comm,ekat::ParameterList());
  auto upgm = std::dynamic_pointer_cast<UserProvidedGridsManager>(gm);
  upgm->set_grid(pg);
  upgm->set_reference_grid(grid_name);

  FieldRepository<Real> repo;
  repo.registration_begins();
  repo.register_field(qv_id,"TRACERS");
  repo.register_field(qc_id,"TRACERS");
  repo.register_field(qr_id,"TRACERS");
  repo.registration_ends(gm);

  auto qv = repo.get_field(qv_id);
  auto qc = repo.get_field(qc_id);
  auto qr = repo.get_field(qr_id);

  auto group = repo.get_field_group("TRACERS",grid_name);
  // The repo should have allocated the group bundled
  REQUIRE (group.m_info->m_bundled);

  const auto& Q_name = group.m_bundle->get_header().get_identifier().name();
  auto Q = repo.get_field(Q_name,grid_name);

  // The bundled field in the group should match the field we get from the repo
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

  upgm->clean_up();
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
