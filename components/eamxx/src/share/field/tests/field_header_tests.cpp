#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_header.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_layout.hpp"
#include "share/field/field_tracking.hpp"

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

  REQUIRE (fl3.dim(COL)==1);
  REQUIRE_THROWS (fl3.dim(CMP));
  REQUIRE (fl3.dim(CMP,false)==3);

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

TEST_CASE("field_header_alias", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  std::vector<FieldTag> tags = {COL, LEV};
  std::vector<int>      dims = {3, 24};

  FieldLayout layout(tags, dims);
  FieldIdentifier fid("f", layout, m/s, "my_grid");

  // Build a committed parent header
  auto parent = std::make_shared<FieldHeader>(fid);
  parent->get_alloc_properties().commit(layout);

  // Set some extra data on the parent
  parent->set_extra_data("key1", 42);

  SECTION ("alias(name)") {
    auto alias = parent->alias("f_alias");

    // Name changed, grid unchanged
    REQUIRE (alias->get_identifier().name()      == "f_alias");
    REQUIRE (alias->get_identifier().get_grid_name() == "my_grid");
    // Layout (tags, dims, names) unchanged
    REQUIRE (alias->get_identifier().get_layout() == parent->get_identifier().get_layout());
    // Shared tracking, alloc props, and extra data
    REQUIRE (alias->get_tracking_ptr()         == parent->get_tracking_ptr());
    REQUIRE (&alias->get_alloc_properties()    == &parent->get_alloc_properties());
    REQUIRE (alias->has_extra_data("key1"));
    REQUIRE (alias->get_extra_data<int>("key1") == 42);
    // Extra data added through alias is visible on parent
    alias->set_extra_data("key2", std::string("hello"));
    REQUIRE (parent->has_extra_data("key2"));
  }

  SECTION ("alias(name, grid_name)") {
    auto alias = parent->alias("f_alias", "other_grid");

    REQUIRE (alias->get_identifier().name()           == "f_alias");
    REQUIRE (alias->get_identifier().get_grid_name()  == "other_grid");
    REQUIRE (alias->get_identifier().get_layout()     == parent->get_identifier().get_layout());
    REQUIRE (alias->get_tracking_ptr()    == parent->get_tracking_ptr());
    REQUIRE (&alias->get_alloc_properties() == &parent->get_alloc_properties());
    REQUIRE (alias->has_extra_data("key1"));
  }

  SECTION ("alias(name, tag_names)") {
    std::map<FieldTag,std::string> tag_names = {{COL,"ncol_d"}, {LEV,"lev_d"}};
    auto alias = parent->alias("f_alias", tag_names);

    REQUIRE (alias->get_identifier().name()          == "f_alias");
    REQUIRE (alias->get_identifier().get_grid_name() == "my_grid");
    // Dim names are renamed in the alias
    REQUIRE (alias->get_identifier().get_layout().name(0) == "ncol_d");
    REQUIRE (alias->get_identifier().get_layout().name(1) == "lev_d");
    // Tags and dims are unchanged
    REQUIRE (alias->get_identifier().get_layout().tags() == tags);
    REQUIRE (alias->get_identifier().get_layout().dims() == dims);
    // Parent layout is NOT altered
    REQUIRE (parent->get_identifier().get_layout().name(0) == "ncol");
    REQUIRE (parent->get_identifier().get_layout().name(1) == "lev");
    // Shared tracking, alloc props, extra data
    REQUIRE (alias->get_tracking_ptr()      == parent->get_tracking_ptr());
    REQUIRE (&alias->get_alloc_properties() == &parent->get_alloc_properties());
    REQUIRE (alias->has_extra_data("key1"));
  }

  SECTION ("alias(name, grid_name, tag_names)") {
    std::map<FieldTag,std::string> tag_names = {{COL,"ncol_d"}};
    auto alias = parent->alias("f_alias", "other_grid", tag_names);

    REQUIRE (alias->get_identifier().name()           == "f_alias");
    REQUIRE (alias->get_identifier().get_grid_name()  == "other_grid");
    REQUIRE (alias->get_identifier().get_layout().name(0) == "ncol_d");
    REQUIRE (alias->get_identifier().get_layout().name(1) == "lev");
    REQUIRE (alias->get_identifier().get_layout().tags() == tags);
    REQUIRE (alias->get_identifier().get_layout().dims() == dims);
    REQUIRE (parent->get_identifier().get_layout().name(0) == "ncol");
    REQUIRE (alias->get_tracking_ptr()      == parent->get_tracking_ptr());
    REQUIRE (&alias->get_alloc_properties() == &parent->get_alloc_properties());
    REQUIRE (alias->has_extra_data("key1"));
  }

  SECTION ("alias of a subfield preserves parent info") {
    // Build a subfield header: slice parent at col=0 along dim 0
    FieldLayout sf_layout({LEV}, {24});
    FieldIdentifier sf_id("sf", sf_layout, m/s, "my_grid");
    auto sf_header = create_subfield_header(sf_id, parent, /*idim=*/0, /*k=*/0, /*dynamic=*/false);

    // Parent info is unaltered
    REQUIRE (parent->get_identifier().get_layout() == layout);

    // Subfield alias (no tag renaming)
    auto sf_alias = sf_header->alias("sf_alias");
    REQUIRE (sf_alias->get_identifier().name() == "sf_alias");
    REQUIRE (sf_alias->get_tracking_ptr()      == sf_header->get_tracking_ptr());
    REQUIRE (&sf_alias->get_alloc_properties() == &sf_header->get_alloc_properties());
    // Parent is preserved in the alias
    REQUIRE (sf_alias->get_parent() != nullptr);
    REQUIRE (sf_alias->get_parent() == sf_header->get_parent());
    // Parent's identifier is NOT altered by the subfield or its alias
    REQUIRE (parent->get_identifier().get_layout() == layout);

    // Subfield alias with tag renaming
    std::map<FieldTag,std::string> tag_names = {{LEV,"lev_d"}};
    auto sf_alias_renamed = sf_header->alias("sf_alias_d", tag_names);
    REQUIRE (sf_alias_renamed->get_identifier().get_layout().name(0) == "lev_d");
    REQUIRE (sf_alias_renamed->get_parent() != nullptr);
    REQUIRE (sf_alias_renamed->get_parent() == sf_header->get_parent());
    // Parent still unaltered
    REQUIRE (parent->get_identifier().get_layout() == layout);
  }
}

} // anonymous namespace
