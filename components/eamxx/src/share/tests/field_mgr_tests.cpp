#include <catch2/catch.hpp>
#include <numeric>

#include "share/manager/field_manager.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_utils.hpp>

namespace scream {

TEST_CASE("field_mgr", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FR  = FieldRequest;
  using SL  = std::list<std::string>;
  using Pack1 = ekat::Pack<Real,8>;
  using Pack2 = ekat::Pack<Real,16>;

  const int ncols1 = 4; const int ncols2 = 7;
  const int nlevs1 = 7; const int nlevs2 = 10;
  const int subview_dim = 1;
  const int subview_slice = 0;

  std::vector<FieldTag> tags1 = {COL,LEV};
  std::vector<FieldTag> tags2 = {COL,CMP,LEV};

  std::vector<FieldTag> bad_tags = {EL,GP,GP};

  std::vector<int> dims1 = {ncols1,nlevs1}; std::vector<int> dims2 = {ncols2,nlevs2};
  std::vector<int> dims3 = {ncols1,10,nlevs1}; std::vector<int> dims4 = {ncols2,10,nlevs2};

  std::vector<int> bad_dims = {2,4,4};

  FID fid1_1("field1", {tags1, dims1},  m/s, "grid1");
  FID fid1_2("field1", {tags1, dims2},  m/s, "grid2");
  FID fid2_1("field2", {tags2, dims3},  m/s, "grid1");
  FID fid2_2("field2", {tags2, dims4},  m/s, "grid2");

  const auto km = 1000*m;
  FID bad1("field1", {tags1, dims1},  m/s, "grid3"); // Bad grid
  FID bad2("field1", {tags1, dims1}, km/s, "grid1"); // Bad units
  FID bad3("field2", {tags1, dims2},  m/s, "grid1"); // Bad layout

  ekat::Comm comm(MPI_COMM_WORLD);
  auto g1 = create_point_grid("grid1",ncols1*comm.size(),nlevs1,comm);
  auto g2 = create_point_grid("grid2",ncols2*comm.size(),nlevs2,comm);
  auto gm = std::make_shared<LibraryGridsManager>(g1, g2);
  FieldManager field_mgr(gm);

  // === Valid registration calls === //
  field_mgr.register_field(FR(fid1_1,Pack1::n));
  field_mgr.register_field(FR{fid1_2,Pack2::n});
  field_mgr.register_field(FR{fid2_1});
  field_mgr.register_field(FR{fid2_1,"group_1"});
  field_mgr.register_field(FR{fid1_2,SL{"group_1", "group_2"}});
  field_mgr.register_field(FR{fid2_2});

  // === Invalid registration calls === //
  REQUIRE_THROWS(field_mgr.register_field(FR{bad1}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));
  REQUIRE_THROWS(field_mgr.register_field(FR{bad2}));

  // Cannot add external fields while registration is happening
  REQUIRE_THROWS(field_mgr.add_field(Field(fid1_1)));

  field_mgr.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(field_mgr.register_field(FR{fid1_1}));

  FID new_fid("new_field", {tags1, dims1},  m/s, "grid1");
  REQUIRE_THROWS (field_mgr.add_field(Field(new_fid))); // Not allocated

  REQUIRE (field_mgr.get_repo("grid1").size()==2);
  REQUIRE (field_mgr.get_repo("grid2").size()==2);

  // Get all fields
  auto f1_1 = field_mgr.get_field(fid1_1);
  auto f1_2 = field_mgr.get_field(fid1_2);
  auto f2_1 = field_mgr.get_field(fid2_1);
  auto f2_2 = field_mgr.get_field(fid2_2);

  // Verify both get_field methods match
  REQUIRE (f1_1 == field_mgr.get_field(fid1_1.name(), fid1_1.get_grid_name()));
  REQUIRE (f1_2 == field_mgr.get_field(fid1_2.name(), fid1_2.get_grid_name()));
  REQUIRE (f2_1 == field_mgr.get_field(fid2_1.name(), fid2_1.get_grid_name()));
  REQUIRE (f2_2 == field_mgr.get_field(fid2_2.name(), fid2_2.get_grid_name()));

  // Try to get invalid fields
  REQUIRE_THROWS(field_mgr.get_field("bad", "grid1"));    // Not in the field_mgr
  REQUIRE_THROWS(field_mgr.get_field(bad1));              // Not in field_mgr
  REQUIRE_THROWS(field_mgr.get_field("field1", "grid3")); // Wrong grid

  // Check that the groups names are in the header. While at it, make sure that case insensitive works fine.
  REQUIRE (ekat::contains(f2_1.get_header().get_tracking().get_groups_names(),"gRouP_1"));
  REQUIRE (ekat::contains(f1_2.get_header().get_tracking().get_groups_names(),"Group_2"));
  REQUIRE (ekat::contains(f1_2.get_header().get_tracking().get_groups_names(),"Group_1"));

  // Check that correct grids requested groups
  REQUIRE (field_mgr.has_group("group_1", "grid1"));
  REQUIRE (field_mgr.has_group("group_1", "grid2"));
  REQUIRE (not field_mgr.has_group("group_2", "grid1"));
  REQUIRE (field_mgr.has_group("group_2", "grid2"));

  // Check that the groups in the field_mgr contain the correct fields
  auto gr1_1 = field_mgr.get_group_info("group_1", "grid1");
  auto gr1_2 = field_mgr.get_group_info("group_1", "grid2");
  auto gr2_2 = field_mgr.get_group_info("group_2", "grid2");
  REQUIRE (gr1_1.m_fields_names.size()==1);
  REQUIRE (gr1_2.m_fields_names.size()==1);
  REQUIRE (gr2_2.m_fields_names.size()==1);
  REQUIRE (ekat::contains(gr1_1.m_fields_names,"Field2"));
  REQUIRE (ekat::contains(gr1_2.m_fields_names,"Field1"));
  REQUIRE (ekat::contains(gr2_2.m_fields_names,"Field1"));

  // Check alloc props for f1 and f2 (which requested pack size > 1)
  auto f1_1_padding = f1_1.get_header().get_alloc_properties().get_padding();
  auto f1_2_padding = f1_2.get_header().get_alloc_properties().get_padding();

  REQUIRE (f1_1_padding==ekat::PackInfo<Pack1::n>::padding(nlevs1));
  REQUIRE (f1_2_padding==ekat::PackInfo<Pack2::n>::padding(nlevs2));

  // Try to subview a field and set the subfield back in the FM
  field_mgr.add_field(f2_1.subfield("field2_1_sf",subview_dim,subview_slice,true));
  REQUIRE (field_mgr.get_repo("grid1").size()==3);

  auto f2_1_sf = field_mgr.get_field("field2_1_sf", "grid1");
  REQUIRE_THROWS (field_mgr.add_field(f2_1_sf)); // Cannot have duplicates
}

TEST_CASE("tracers_group", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FR  = FieldRequest;

  const int ncols1 = 4;
  const int ncols2 = 3;
  const int nlevs = 7;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims1 = {ncols1,nlevs};
  std::vector<int> dims2 = {ncols2,nlevs};

  const auto nondim = Units::nondimensional();

  const std::string gn1 = "grid1";
  const std::string gn2 = "grid2";

  FieldIdentifier qv_id("qv", {tags, dims1}, nondim, gn1);
  FieldIdentifier a_id("a", {tags, dims1}, nondim, gn1);
  FieldIdentifier b_id("b", {tags, dims2}, nondim, gn2);
  FieldIdentifier c_id("c", {tags, dims1}, nondim, gn1);

  ekat::Comm comm(MPI_COMM_WORLD);
  auto g1 = create_point_grid(gn1,ncols1*comm.size(),nlevs,comm);
  auto g2 = create_point_grid(gn2,ncols2*comm.size(),nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>(g1, g2);
  FieldManager field_mgr(gm);

  using los = std::list<std::string>;
  field_mgr.register_field(FR{qv_id,"tracers"});
  field_mgr.register_field(FR{a_id,"tracers"});
  field_mgr.register_field(FR{b_id,los{"tracers", "subtracers"}});
  field_mgr.register_field(FR{c_id,los{"tracers", "subtracers"}});

  field_mgr.register_group(GroupRequest("tracers",gn1,MonolithicAlloc::Required));
  field_mgr.register_group(GroupRequest("tracers",gn2,MonolithicAlloc::Required));
  field_mgr.register_group(GroupRequest("subtracers",gn1,MonolithicAlloc::Required));
  //field_mgr.register_group(GroupRequest("subtracers",gn2,MonolithicAlloc::Required));

  field_mgr.registration_ends();

  auto T1 = field_mgr.get_field("tracers", gn1);
  auto qv1 = field_mgr.get_field("qv", gn1);
  auto a1 = field_mgr.get_field("a", gn1);
  auto b1 = field_mgr.get_field("b", gn1);
  auto c1 = field_mgr.get_field("c", gn1);
  auto T2 = field_mgr.get_field("tracers", gn2);
  auto qv2 = field_mgr.get_field("qv", gn2);
  auto a2 = field_mgr.get_field("a", gn2);
  auto b2 = field_mgr.get_field("b", gn2);
  auto c2 = field_mgr.get_field("c", gn2);

  // No subtracer field should exist since it is not
  // a parent of any field.
  REQUIRE_THROWS (field_mgr.get_field("subtracers", gn1));
  REQUIRE_THROWS (field_mgr.get_field("subtracers", gn2));

  // The field_mgr should have allocated the group as a monolith
  auto tracers1 = field_mgr.get_field_group("tracers", gn1);
  auto tracers2 = field_mgr.get_field_group("tracers", gn2);
  auto subtracers = field_mgr.get_field_group("subtracers", gn1);
  REQUIRE (tracers1.m_info->m_monolithic_allocation);
  REQUIRE (tracers2.m_info->m_monolithic_allocation);
  REQUIRE (subtracers.m_info->m_monolithic_allocation);

  // The monolithic field in the tracers group should match the field we get from the field_mgr
  REQUIRE (T1.is_aliasing(*tracers1.m_monolithic_field));
  REQUIRE (T2.is_aliasing(*tracers2.m_monolithic_field));

  // Require that the parent of each field is the "tracers" group
  auto qv1_p = qv1.get_header().get_parent();
  auto a1_p = a1.get_header().get_parent();
  auto b1_p = b1.get_header().get_parent();
  auto c1_p = c1.get_header().get_parent();
  auto qv2_p = qv2.get_header().get_parent();
  auto a2_p = a2.get_header().get_parent();
  auto b2_p = b2.get_header().get_parent();
  auto c2_p = c2.get_header().get_parent();

  REQUIRE ((qv1_p!=nullptr && qv1_p.get()==&T1.get_header()));
  REQUIRE ((a1_p!=nullptr && a1_p.get()==&T1.get_header()));
  REQUIRE ((b1_p!=nullptr && b1_p.get()==&T1.get_header()));
  REQUIRE ((c1_p!=nullptr && c1_p.get()==&T1.get_header()));
  REQUIRE ((qv2_p!=nullptr && qv2_p.get()==&T2.get_header()));
  REQUIRE ((a2_p!=nullptr && a2_p.get()==&T2.get_header()));
  REQUIRE ((b2_p!=nullptr && b2_p.get()==&T2.get_header()));
  REQUIRE ((c2_p!=nullptr && c2_p.get()==&T2.get_header()));

  // Require subtracers monolith is subfield of tracers
  REQUIRE ((
    subtracers.m_monolithic_field->get_header().get_parent() != nullptr &&
    subtracers.m_monolithic_field->get_header().get_parent().get() == &T1.get_header()));

  const auto idx_qv1 = tracers1.m_info->m_subview_idx.at("qv");
  const auto idx_a1 = tracers1.m_info->m_subview_idx.at("a");
  const auto idx_b1 = tracers1.m_info->m_subview_idx.at("b");
  const auto idx_c1 = tracers1.m_info->m_subview_idx.at("c");
  const auto idx_qv2 = tracers2.m_info->m_subview_idx.at("qv");
  const auto idx_a2 = tracers2.m_info->m_subview_idx.at("a");
  const auto idx_b2 = tracers2.m_info->m_subview_idx.at("b");
  const auto idx_c2 = tracers2.m_info->m_subview_idx.at("c");

  const auto sub_idx_b1 = subtracers.m_info->m_subview_idx.at("b");
  const auto sub_idx_c1 = subtracers.m_info->m_subview_idx.at("c");

  // Subview indices of all groups should be identical over requested grids
  REQUIRE (idx_qv1 == idx_qv2);
  REQUIRE (idx_a1 == idx_a2);
  REQUIRE (idx_b1 == idx_b2);
  REQUIRE (idx_c1 == idx_c2);

  // Indices should be in the interval [0,4)
  REQUIRE ((0<=idx_qv1 and idx_qv1<=3));
  REQUIRE ((0<=idx_a1 and idx_a1<=3));
  REQUIRE ((0<=idx_b1 and idx_b1<=3));
  REQUIRE ((0<=idx_c1 and idx_c1<=3));

  // Indices of sub tracers should be between [0,2)
  REQUIRE ((0<=sub_idx_b1 and sub_idx_b1<=1));
  REQUIRE ((0<=sub_idx_c1 and sub_idx_c1<=1));

  // Now fill tracer field with random values
  auto engine = setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0,1.0);

  randomize(T1,engine,pdf);
  randomize(T2,engine,pdf);

  // Check that the same values are in all individual tracers
  T1.sync_to_host();
  auto T1_h = T1.get_view<Real***,Host>();
  auto qv1_h = qv1.get_view<Real**,Host>();
  auto a1_h = a1.get_view<Real**,Host>();
  auto b1_h = b1.get_view<Real**,Host>();
  auto c1_h = c1.get_view<Real**,Host>();

  for (int icol=0; icol<ncols1; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      REQUIRE (T1_h(icol,idx_qv1,ilev)==qv1_h(icol,ilev));
      REQUIRE (T1_h(icol,idx_a1,ilev)==a1_h(icol,ilev));
      REQUIRE (T1_h(icol,idx_b1,ilev)==b1_h(icol,ilev));
      REQUIRE (T1_h(icol,idx_c1,ilev)==c1_h(icol,ilev));
    }
  }

  T2.sync_to_host();
  auto T2_h = T2.get_view<Real***,Host>();
  auto qv2_h = qv2.get_view<Real**,Host>();
  auto a2_h = a2.get_view<Real**,Host>();
  auto b2_h = b2.get_view<Real**,Host>();
  auto c2_h = c2.get_view<Real**,Host>();

  for (int icol=0; icol<ncols2; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      REQUIRE (T2_h(icol,idx_qv2,ilev)==qv2_h(icol,ilev));
      REQUIRE (T2_h(icol,idx_a2,ilev)==a2_h(icol,ilev));
      REQUIRE (T2_h(icol,idx_b2,ilev)==b2_h(icol,ilev));
      REQUIRE (T2_h(icol,idx_c2,ilev)==c2_h(icol,ilev));
    }
  }

  // Check that changing sub tracers change tracers group and individual tracers
  T1.deep_copy(0.0);
  auto& sub_T1 = subtracers.m_monolithic_field;
  randomize(*sub_T1,engine,pdf);
  sub_T1->sync_to_host();
  auto sub_T1_h = sub_T1->get_strided_view<Real***,Host>();

  for (int icol=0; icol<ncols1; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      REQUIRE (sub_T1_h(icol,sub_idx_b1,ilev)==b1_h(icol,ilev));
      REQUIRE (sub_T1_h(icol,sub_idx_c1,ilev)==c1_h(icol,ilev));
    }
  }

  // Check that the field ptrs stored in the group are the same as the fields
  auto qv1_ptr = tracers1.m_individual_fields.at("qv");
  auto a1_ptr = tracers1.m_individual_fields.at("a");
  auto b1_ptr = tracers1.m_individual_fields.at("b");
  auto c1_ptr = tracers1.m_individual_fields.at("c");

  REQUIRE (qv1_ptr->is_aliasing(qv1));
  REQUIRE (a1_ptr->is_aliasing(a1));
  REQUIRE (b1_ptr->is_aliasing(b1));
  REQUIRE (c1_ptr->is_aliasing(c1));

  auto qv2_ptr = tracers2.m_individual_fields.at("qv");
  auto a2_ptr = tracers2.m_individual_fields.at("a");
  auto b2_ptr = tracers2.m_individual_fields.at("b");
  auto c2_ptr = tracers2.m_individual_fields.at("c");

  REQUIRE (qv2_ptr->is_aliasing(qv2));
  REQUIRE (a2_ptr->is_aliasing(a2));
  REQUIRE (b2_ptr->is_aliasing(b2));
  REQUIRE (c2_ptr->is_aliasing(c2));

  b1_ptr = subtracers.m_individual_fields.at("b");
  c1_ptr = subtracers.m_individual_fields.at("c");
  REQUIRE (b1_ptr->is_aliasing(b1));
  REQUIRE (c1_ptr->is_aliasing(c1));
}

} // anonymous scream
