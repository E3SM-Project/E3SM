#include "catch2/catch.hpp"
#include "physics/shoc/shoc_f90.hpp"
#include "physics/shoc/shoc_ic_cases.hpp"

namespace {

TEST_CASE("FortranData", "shoc") {
  int val = scream::shoc::test_FortranData();
  REQUIRE(val == 0);
}

TEST_CASE("FortranDataIterator", "shoc") {
  using scream::shoc::ic::Factory;
  const auto d = Factory::create();
  scream::shoc::FortranDataIterator fdi(d);
  REQUIRE(fdi.nfield() == 43);
  const auto& f = fdi.getfield(0);
  REQUIRE(f.dim == 2);
  REQUIRE(f.extent[0] == d->shcol);
  REQUIRE(f.extent[1] == 1);
  REQUIRE(f.extent[2] == 1);
  REQUIRE(f.data == d->host_dx.data());
  REQUIRE(f.size == d->shcol);
}

TEST_CASE("shoc_init_f", "shoc") {
  int nerr = scream::shoc::test_shoc_init(true);
  REQUIRE(nerr == 0);
}

TEST_CASE("shoc_ic_f", "shoc") {
  int nerr = scream::shoc::test_shoc_ic(true);
  REQUIRE(nerr == 0);
}

TEST_CASE("shoc_init_c", "shoc") {
  int nerr = scream::shoc::test_shoc_init(false);
  REQUIRE(nerr == 0);
}

TEST_CASE("shoc_ic_c", "shoc") {
  int nerr = scream::shoc::test_shoc_ic(false);
  REQUIRE(nerr == 0);
}


} // empty namespace
