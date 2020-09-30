#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"

#include "ekat/ekat_pack.hpp"

namespace {

TEST_CASE("field_identifier", "") {
  using namespace scream;
  using namespace ekat::units;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::VerticalLevel};

  FieldIdentifier fid1 ("field_1", tags1, m/s);
  FieldIdentifier fid2 ("field_1", tags1, m/s);
  FieldIdentifier fid3 ("field_1", tags2, m/s);
  FieldIdentifier fid4 ("field_2", tags2, m/s);

  REQUIRE (fid1==fid2);
  REQUIRE (fid2!=fid3);
  REQUIRE (fid3!=fid4);

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};

  // Should not be able to set negative dimensions, or a non existing extent
  REQUIRE_THROWS(fid1.set_dimension(0,-1));
  REQUIRE_THROWS(fid1.set_dimension(-1,1));
  REQUIRE_THROWS(fid1.set_dimension(3,1));

  // Should not be able to set a dimensions vector of wrong rank
  REQUIRE_THROWS(fid1.set_dimensions({1,2}));

  fid1.set_dimensions(dims1);
  fid2.set_dimensions(dims2);

  // Should not be able to reset the dimensions once they are set
  REQUIRE_THROWS(fid2.set_dimensions(dims2));

  REQUIRE (fid1!=fid2);
}

TEST_CASE("field_tracking", "") {
  using namespace scream;

  FieldTracking track("track");
  util::TimeStamp time1(0,0,0,10.0);
  util::TimeStamp time2(0,0,0,20.0);
  REQUIRE_NOTHROW (track.update_time_stamp(time2));

  REQUIRE_THROWS  (track.update_time_stamp(time1));
}
  template<typename ST, int N>
  using Pack = ekat::Pack<ST,N>;

TEST_CASE("field", "") {

  using namespace scream;
  using namespace ekat::units;

  using Device = DefaultDevice;

  std::vector<FieldTag> tags = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::VerticalLevel};
  std::vector<int> dims = {2, 3, 12};

  FieldIdentifier fid ("field_1", tags, m/s);
  fid.set_dimensions(dims);

  // Check copy constructor
  SECTION ("copy ctor") {
    Field<Real,Device> f1 (fid);

    f1.allocate_view();

    Field<const Real,Device> f2 = f1;
    REQUIRE(f2.get_header_ptr()==f1.get_header_ptr());
    REQUIRE(f2.get_view()==f1.get_view());
    REQUIRE(f2.is_allocated());
  }

  // Check if we can extract a reshaped view
  SECTION ("reshape simple") {
    Field<Real,Device> f1 (fid);

    // Should not be able to reshape before allocating
    REQUIRE_THROWS(f1.get_reshaped_view<Real*>());

    f1.allocate_view();

    // Should not be able to reshape to this data type
    REQUIRE_THROWS(f1.get_reshaped_view<Pack<Real,8>>());

    auto v1d = f1.get_view();
    auto v3d = f1.get_reshaped_view<Real[2][3][12]>();
    REQUIRE(v3d.size()==v1d.size());
  }

  // Check if we can request multiple value types
  SECTION ("reshape multiple value types") {
    Field<Real,Device> f1 (fid);
    f1.get_header().get_alloc_properties().request_value_type_allocation<Pack<Real,8>>();
    f1.allocate_view();

    auto v1d = f1.get_view();
    auto v3d_1 = f1.get_reshaped_view<Pack<Real,8>***>();
    auto v3d_2 = f1.get_reshaped_view<Pack<Real,4>***>();
    auto v3d_3 = f1.get_reshaped_view<Real***>();
    auto v3d_4 = f1.get_reshaped_view<Real[2][3][16]>();

    // The memory spans should be identical
    REQUIRE (v3d_1.impl_map().memory_span()==v3d_2.impl_map().memory_span());
    REQUIRE (v3d_1.impl_map().memory_span()==v3d_3.impl_map().memory_span());
    REQUIRE (v3d_1.impl_map().memory_span()==v3d_4.impl_map().memory_span());

    // Sizes differ, since they are in terms of the stored value type.
    // Each Pack<Real,8> corresponds to two Pack<Real,4>, which corresponds to 4 Real's.
    REQUIRE(2*v3d_1.size()==v3d_2.size());
    REQUIRE(8*v3d_1.size()==v3d_3.size());
    REQUIRE(8*v3d_1.size()==v3d_4.size());

    // Trying to reshape into something that the allocation cannot accommodate should throw
    REQUIRE_THROWS (f1.get_reshaped_view<Pack<Real,32>***>());
  }
}

TEST_CASE("field_repo", "") {
  using namespace scream;
  using namespace ekat::units;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Column};

  const auto km = 1000*m;

  FieldIdentifier fid1("field_1", tags1, m/s);
  FieldIdentifier fid2("field_2", tags1, m/s);
  FieldIdentifier fid3("field_2", tags2, m/s);
  FieldIdentifier fid4("field_2", tags2, km/s);

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};
  std::vector<int> dims3 = {13};

  fid1.set_dimensions(dims1);
  fid2.set_dimensions(dims2);
  fid3.set_dimensions(dims3);
  fid4.set_dimensions(dims3);

  FieldRepository<Real,DefaultDevice>  repo;

  // Should not be able to register fields yet
  REQUIRE_THROWS(repo.register_field(fid1,"group_1"));

  repo.registration_begins();
  repo.register_field(fid1,"group_1");
  repo.register_field(fid2,"group_2");
  repo.register_field(fid3,"group_2");
  // Should not be able to register fields to the 'state' group (it's reserved)
  REQUIRE_THROWS(repo.register_field(fid2,"state"));
  // Should not be able to register the same field name with two different units
  REQUIRE_THROWS(repo.register_field(fid4));
  repo.registration_ends();

  // Should not be able to register fields anymore
  REQUIRE_THROWS(repo.register_field(fid1,"group_1"));

  // Check registration is indeed closed
  REQUIRE (repo.repository_state()==RepoState::Closed);
  REQUIRE (repo.size()==2);
  REQUIRE (repo.internal_size()==3);

  auto f1 = repo.get_field(fid1);
  auto f2 = repo.get_field(fid2);

  // Check the two fields identifiers are indeed different
  REQUIRE (f1.get_header().get_identifier()!=f2.get_header().get_identifier());

  // Check that the groups names are in the header. While at it, make sure that case insensitive works fine.
  REQUIRE (ekat::contains(f1.get_header().get_tracking().get_groups_names(),"gRouP_1"));
  REQUIRE (ekat::contains(f2.get_header().get_tracking().get_groups_names(),"Group_2"));

  // Check that the groups in the repo contain the correct fields
  REQUIRE (repo.get_field_groups().count("GROUP_1")==1);
  REQUIRE (repo.get_field_groups().count("GRoup_2")==1);
  REQUIRE (ekat::contains(repo.get_field_groups().at("group_1"),"Field_1"));
  REQUIRE (ekat::contains(repo.get_field_groups().at("group_2"),"Field_2"));
}

} // anonymous namespace
