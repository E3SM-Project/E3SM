#include "catch2/catch.hpp"

#include <share/field/field_identifier.hpp>
#include <share/field/field_header.hpp>
#include <share/field/field.hpp>
#include <share/field/field_repository.hpp>
#include <share/field/field_utils.hpp>
#include <share/scream_pack.hpp>

namespace {

TEST_CASE("field_identifier", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::Level};

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

  std::vector<FieldTag> tags = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::Level};
  std::vector<int> dims1 = {2, 3, 4};

  FieldIdentifier fid1 ("field_1", tags);
  fid1.set_dimensions(dims1);

  // The following 'REQUIRE' tests are not the real tests. The real
  // tests are whether or not the instructions to create fields succeed.
  // If they fail, they'll raise MPI_Abort errors.

  Field<Real*,HostMemSpace,MemoryManaged> f1 (fid1);
  REQUIRE(static_cast<bool>(f1.get_header_ptr()));
  f1.allocate_view();

  Field<const Real*,HostMemSpace,MemoryManaged> f2 = f1;
  REQUIRE(static_cast<bool>(f2.get_header_ptr()));

  auto f3 = reinterpret_field<Real[2][3][4]>(f1);
  REQUIRE(static_cast<bool>(f3.get_header_ptr()));

  // Check if we can copy f3 to a field with a different scalar type (like a Pack)
  FieldIdentifier fid2 ("field_1",tags);
  std::vector<int> dims2 = {2, 3, 1};
  fid2.set_dimensions(dims2);

  auto f3_pack = reinterpret_field<pack::Pack<Real,4>[2][3][1]>(f3);
  const int extent3_f3_pack = f3_pack.get_view().extent_int(2);
  const int extent3_f3      = f3.get_view().extent_int(2);
  REQUIRE(extent3_f3==4);
  REQUIRE(extent3_f3_pack==1);
}

TEST_CASE("field_repo", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::Level};

  FieldIdentifier fid1("field_1", tags1);
  FieldIdentifier fid2("field_2", tags1);

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};

  fid1.set_dimensions(dims1);
  fid2.set_dimensions(dims2);

  FieldRepository<Real,ExecMemSpace>  repo_dev;
  repo_dev.register_field(fid1);
  repo_dev.register_field(fid2);
  repo_dev.registration_complete();

  auto f1 = repo_dev.get_field(fid1);
  auto f2 = repo_dev.get_field(fid2);

  // Check the two fields identifiers are indeed different
  REQUIRE (f1.get_header().get_identifier()!=f2.get_header().get_identifier());

  // Check registration is indeed closed
  REQUIRE (!repo_dev.is_registration_open());
}

} // anonymous namespace
