#include "catch2/catch.hpp"
#include "p3_f90.hpp"
#include "p3_main_wrap.hpp"
#include "p3_ic_cases.hpp"

namespace {

TEST_CASE("FortranData", "p3") {
  int val = scream::p3::test_FortranData();
  REQUIRE(val == 0);
}

TEST_CASE("FortranDataIterator", "p3") {
  using scream::p3::ic::Factory;
  const auto d = Factory::create(Factory::mixed);
  scream::p3::FortranDataIterator fdi(d);
  REQUIRE(fdi.nfield() == 34);
  const auto& f = fdi.getfield(0);
  REQUIRE(f.dim == 2);
  REQUIRE(f.extent[0] == 1);
  REQUIRE(f.extent[1] == 72);
  REQUIRE(f.extent[2] == 1);
  REQUIRE(f.data == d->qv.data());
  REQUIRE(f.size == 72);
}

TEST_CASE("p3_init", "p3") {
  int nerr = scream::p3::test_p3_init();
  REQUIRE(nerr == 0);
}

TEST_CASE("p3_ic_f", "p3") {
  int nerr = scream::p3::test_p3_ic(true);
  REQUIRE(nerr == 0);
}

TEST_CASE("p3_ic_c", "p3") {
  int nerr = scream::p3::test_p3_ic(false);
  REQUIRE(nerr == 0);
}

} // empty namespace
