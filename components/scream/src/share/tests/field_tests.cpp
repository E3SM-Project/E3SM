#include "catch2/catch.hpp"

#include "share/field_header.hpp"
#include "share/field.hpp"
#include "share/field_repository.hpp"

namespace {

TEST_CASE("field_header", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::Level};
  std::vector<FieldTag> tags3 = {FieldTag::Column, FieldTag::Level};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};
  std::vector<int> dims3 = {100, 72};

  FieldHeader fh1 ("field_1", tags1, {2, 3, 4});
  FieldHeader fh2 ("field_1", tags1);
  fh2.set_dimensions(dims1);
  FieldHeader fh3 ("field_1", tags2, dims1);
  FieldHeader fh4 ("field_2", tags2, dims1);
  FieldHeader fh5 ("field_2", tags2);
  fh5.set_dimensions({2,3,3});

  REQUIRE (fh1==fh2);
  REQUIRE (fh2!=fh3);
  REQUIRE (fh3!=fh4);
  REQUIRE (fh4!=fh5);
}

TEST_CASE("field_repo", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::Level};
  std::vector<FieldTag> tags3 = {FieldTag::Column, FieldTag::Level};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};
  std::vector<int> dims3 = {100, 72};

  std::shared_ptr<FieldHeader> fh1 = std::make_shared<FieldHeader>("field_1", tags1, dims1);
  std::shared_ptr<FieldHeader> fh2 = std::make_shared<FieldHeader>("field_1", tags1, dims2);
  std::shared_ptr<FieldHeader> fh3 = std::make_shared<FieldHeader>("field_2", tags1, dims1);
  std::shared_ptr<FieldHeader> fh4 = std::make_shared<FieldHeader>("field_2", tags2, dims1);

  FieldRepository<ExecMemSpace>  repo_dev;
  auto f1 = repo_dev.register_field(fh1);
  auto f2 = repo_dev.register_field(fh2);
  repo_dev.registration_complete();

  // Check the two fields headers are different
  REQUIRE (f1.header()!=f2.header());

  // Reshape a field, and check it went well;
  Field<Real[2][3][4],ExecMemSpace,false> f1_shaped(f1);
  REQUIRE (f1_shaped.header()==f1.header());
  REQUIRE (decltype(f1)::view_type::Rank==1);
  REQUIRE (decltype(f1_shaped)::view_type::Rank==3);
  REQUIRE (f1.view().size()==f1_shaped.view().size());
}

} // anonymous namespace
