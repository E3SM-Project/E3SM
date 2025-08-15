#include <catch2/catch.hpp>

#include "share/util/eamxx_combine_ops.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/eamxx_types.hpp"

#include <ekat_pack.hpp>

#include <iomanip>

namespace {

TEST_CASE ("combine_ops") {
  using namespace scream;
  using pack_type = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  constexpr auto Replace   = CombineMode::Replace;
  constexpr auto Update    = CombineMode::Update;
  constexpr auto Multiply  = CombineMode::Multiply;
  constexpr auto Divide    = CombineMode::Divide;
  constexpr auto Max       = CombineMode::Max;
  constexpr auto Min       = CombineMode::Min;
  constexpr auto fv_val    = constants::fill_value<Real>;

  const pack_type two (2.0);
  const pack_type four (4.0);
  const pack_type six (6.0);
  const pack_type ten (10.0);
  const pack_type fv (constants::fill_value<Real>);
  const pack_type fv_times_ten (constants::fill_value<Real> * 10);

  pack_type x;

  x = two;
  combine<Replace>(two,x,1,0);
  REQUIRE ( (x==two).all() );

  combine<Update>(two,x,2.0,1.0);
  REQUIRE ( (x==six).all() );
  fill_aware_combine<Update>(fv,x,fv_val,2.0,1.0);
  if (not (x==six).all() ) {
    std::cout << "x: " << x << "\n";
    std::cout << " x[0]: " << std::setprecision(18) << x[0] << "\n";
    std::cout << "fv[0]: " << std::setprecision(18) << fv[0] << "\n";
  }
  REQUIRE ( (x==six).all() );

  x = two;
  combine<Multiply>(two,x,1,1);
  REQUIRE ( (x==four).all() );
  fill_aware_combine<Multiply>(fv,x,fv_val,1,1);
  REQUIRE ( (x==four).all() );

  x = four;
  combine<Divide>(two,x,1,1);
  REQUIRE ( (x==two).all() );
  fill_aware_combine<Divide>(fv,x,fv_val,1,1);
  REQUIRE ( (x==two).all() );

  x = two;
  combine<Max>(four,x,1,1);
  REQUIRE ( (x==four).all() );
  fill_aware_combine<Max>(fv,x,fv_val,1,1);
  REQUIRE ( (x==four).all() );

  x = four;
  combine<Min>(two,x,1,1);
  REQUIRE ( (x==two).all() );
  x = fv_times_ten;
  fill_aware_combine<Min>(fv,x,fv_val,1,1);
  REQUIRE ( (x==fv_times_ten).all() );
}

} // anonymous namespace
