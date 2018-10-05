#include "catch2/catch.hpp"

#include "share/field_header.hpp"

namespace {

TEST_CASE("field_header", "") {
  using namespace scream;

  std::vector<FieldTag> tags1 = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
  std::vector<FieldTag> tags2 = {FieldTag::Element, FieldTag::Component, FieldTag::Level};
  std::vector<FieldTag> tags3 = {FieldTag::Column, FieldTag::Level};

  std::vector<int> dims1 = {2, 3, 4};
  std::vector<int> dims2 = {2, 3, 3};
  std::vector<int> dims3 = {100, 72};

  FieldHeader fh1 ("field_1", {2, 3, 4}, tags1);
  FieldHeader fh2 ("field_1", tags1);
  fh2.set_dimensions(dims1);
  FieldHeader fh3 ("field_1", dims1, tags2);
  FieldHeader fh4 ("field_2", dims1, tags2);
  FieldHeader fh5 ("field_2", tags2);
  fh5.set_dimensions({2,3,3});

  REQUIRE (fh1==fh2);
  REQUIRE (fh2!=fh3);
  REQUIRE (fh3!=fh4);
  REQUIRE (fh4!=fh5);
}

} // anonymous namespace

