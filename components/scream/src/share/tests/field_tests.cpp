#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"
#include "share/scream_pack.hpp"

namespace {

TEST_CASE("field_identifier", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::VerticalLevel};

  FieldIdentifier fid1 ("field_1", tags1);
  FieldIdentifier fid2 ("field_1", tags1);
  FieldIdentifier fid3 ("field_1", tags2);
  FieldIdentifier fid4 ("field_2", tags2);

  REQUIRE (fid1==fid2);
  REQUIRE (fid2!=fid3);
  REQUIRE (fid3!=fid4);

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};

  fid1.set_dimensions(dims1);
  fid2.set_dimensions(dims2);

  REQUIRE (fid1!=fid2);
}

TEST_CASE("field", "") {
  using namespace scream;
  using namespace scream::pack;

  using Device = DefaultDevice;

  std::vector<FieldTag> tags = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::VerticalLevel};
  std::vector<int> dims = {2, 3, 12};

  FieldIdentifier fid ("field_1", tags);
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
    f1.allocate_view();

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
  }
}

TEST_CASE("field_repo", "") {
  using namespace scream;

  using Device = DefaultDevice;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::VerticalLevel};

  FieldIdentifier fid1("field_1", tags1);
  FieldIdentifier fid2("field_2", tags1);

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};

  fid1.set_dimensions(dims1);
  fid2.set_dimensions(dims2);

  FieldRepository<Real,Device>  repo_dev;

  repo_dev.registration_begins();
  repo_dev.register_field(fid1);
  repo_dev.register_field(fid2);
  repo_dev.registration_ends();

  // Check registration is indeed closed
  REQUIRE (repo_dev.repository_state()==RepoState::Closed);
  REQUIRE (repo_dev.size()==2);

  auto f1 = repo_dev.get_field(fid1);
  auto f2 = repo_dev.get_field(fid2);

  // Check the two fields identifiers are indeed different
  REQUIRE (f1.get_header().get_identifier()!=f2.get_header().get_identifier());
}

} // anonymous namespace
