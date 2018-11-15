#include "catch2/catch.hpp"

#include <share/field/field_identifier.hpp>
#include <share/field/field_header.hpp>
#include <share/field/field.hpp>
#include <share/field/field_repository.hpp>

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

  FieldRepository<ExecMemSpace>  repo_dev;
  auto f1 = repo_dev.register_field(fid1);
  auto f2 = repo_dev.register_field(fid2);
  repo_dev.registration_complete();

  // Check the two fields identifiers are indeed different
  REQUIRE (f1.get_header().get_identifier()!=f2.get_header().get_identifier());

  // Check registration is indeed closed
  REQUIRE (!repo_dev.is_registration_open());
}

} // anonymous namespace
