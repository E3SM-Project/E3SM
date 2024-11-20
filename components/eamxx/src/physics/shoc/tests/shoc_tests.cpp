#include "catch2/catch.hpp"

#include "shoc_main_wrap.hpp"
#include "shoc_ic_cases.hpp"

#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("FortranData", "shoc") {
  int val = scream::shoc::test_FortranData();
  REQUIRE(val == 0);
}

TEST_CASE("FortranDataIterator", "shoc") {
  using scream::shoc::ic::Factory;
  const auto d = Factory::create(Factory::standard);
  scream::shoc::FortranDataIterator fdi(d);
  REQUIRE(fdi.nfield() == 44);
  const auto& f = fdi.getfield(0);
  REQUIRE(f.dim == 2);
  REQUIRE(f.extent[0] == d->shcol);
  REQUIRE(f.extent[1] == 1);
  REQUIRE(f.extent[2] == 1);
  REQUIRE(f.data == d->host_dx.data());
  REQUIRE(static_cast<int>(f.size) == d->shcol);
}

TEST_CASE("shoc_ic_c", "shoc") {
  int nerr = scream::shoc::test_shoc_ic();
  REQUIRE(nerr == 0);
}

} // anonymous namespace

