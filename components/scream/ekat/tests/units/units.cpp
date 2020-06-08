#include <catch2/catch.hpp>

#include "ekat/util/scream_units.hpp"

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
#if defined(SCREAM_CONSTEXPR_ASSERT) && !defined(NDEBUG)
    const RationalConstant zero(0);
    REQUIRE_THROWS(one/zero);
    REQUIRE_THROWS(pow(zero,zero));
    REQUIRE_THROWS(pow(-one,half));
#endif
  }

  SECTION ("scaling_factor") {
    constexpr RationalConstant one = RationalConstant::one();
    constexpr RationalConstant third = 1/(3*one);
    constexpr RationalConstant four_thirds = 4*third;
    constexpr RationalConstant three_halves = 3/(2*one);
    constexpr RationalConstant two = 2*one;
    constexpr ScalingFactor root2 (2,RationalConstant{1,2});

    // I need a runtime number for runtime checks
    const RationalConstant three(3);

    // Verify operations
    REQUIRE( root2*root2 == two );
    REQUIRE( sqrt(two) == root2 );
    REQUIRE( pow(three_halves,2)*pow(four_thirds,3) == 16*one/3);

    // Verify printing
    REQUIRE(to_string(root2)=="2^1/2");
    REQUIRE(to_string(root2,Format::Float)=="2^0.5");

#if defined(SCREAM_CONSTEXPR_ASSERT) && !defined(NDEBUG)
    REQUIRE_THROWS(sqrt(-three));
#endif
  }

  SECTION ("units") {
    constexpr RationalConstant one = RationalConstant::one();
    const auto km = kilo*m;
    const auto kPa = kilo*Pa;

    Units nondim (ScalingFactor(1));
    Units milliJ = milli*N*m;
    Units mix_ratio = kg/kg;
    mix_ratio.set_string("kg/kg");

    // Verify operations
    REQUIRE (milliJ == kPa*pow(m,3)/mega);
    REQUIRE (m/s*day/km == Units(one*86400/1000));
    REQUIRE (pow(sqrt(m),2)==m);

    // Verify printing
    REQUIRE (to_string(nondim)=="1");
    REQUIRE (to_string(milliJ.set_exp_format(Format::Rat))=="0.001 m^2 s^-2 kg");

    // Verify changing the string works and does not affect the to_string function
    REQUIRE (mix_ratio==nondim);
    REQUIRE (to_string(mix_ratio)=="1");
    REQUIRE (mix_ratio.get_string()=="kg/kg");
  }
}
