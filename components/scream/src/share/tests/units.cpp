#include <catch2/catch.hpp>

#include "share/util/units.hpp"

#include <iostream>

TEST_CASE("units_framework", "") {
  using namespace scream;
  using namespace scream::units;

  SECTION ("rational_constant") {
    constexpr RationalConstant quarter(1,4);
    constexpr RationalConstant half(1,2);
    constexpr RationalConstant one(1);
    constexpr RationalConstant two(2,1);

    // Verify operations
    REQUIRE(half*two == one);
    REQUIRE(half+half == one);
    REQUIRE(one/two == half);
    REQUIRE(half*half == quarter);
    REQUIRE(pow(half,2)==quarter);
    REQUIRE_THROWS(RationalConstant(1,0));
    REQUIRE_THROWS(RationalConstant(0,0));
    REQUIRE_THROWS(pow(0*one,0));
  }

  SECTION ("rational_constant") {
    constexpr RationalConstant one = RationalConstant::one();
    constexpr RationalConstant third = 1/(3*one);
    constexpr RationalConstant four_thirds = 4*third;
    constexpr RationalConstant three_halves = 3/(2*one);
    constexpr RationalConstant two = 2*one;
    constexpr ScalingFactor root2 (2,{1,2});

    // Verify operations
    REQUIRE( root2*root2 == two );
    REQUIRE( sqrt(two) == root2 );
    REQUIRE( pow(three_halves,2)*pow(four_thirds,3) == 16*one/3);

    // Verify printing
    REQUIRE(to_string(root2)=="2^1/2");
    REQUIRE(to_string(root2,Format::Float)=="2^0.5");
  }

  SECTION ("units") {
    constexpr RationalConstant one = RationalConstant::one();
    constexpr auto km = kilo*m;
    constexpr auto kPa = kilo*Pa;

    Units nondim (ScalingFactor(1));
    Units milliJ = milli*N*m;

    // Verify operations
    REQUIRE (milliJ == kPa*pow(m,3)/mega);
    REQUIRE (m/s*day/km == Units(one*86400/1000));
    REQUIRE (pow(sqrt(m),2)==m);

    // Verify printing
    REQUIRE (to_string(nondim)=="1");
    REQUIRE (to_string(milliJ.set_exp_format(Format::Rat))=="0.001 m^2 s^-2 kg");
  }
}
