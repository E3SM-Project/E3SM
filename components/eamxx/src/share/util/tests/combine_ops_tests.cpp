#include <catch2/catch.hpp>

#include "share/util/eamxx_combine_ops.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/core/eamxx_types.hpp"

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

  combine<Update>(two,x,sp(2.0),sp(1.0));
  REQUIRE ( (x==six).all() );
  combine_fill_aware<Update>(fv,x,sp(2.0),sp(1.0));
  REQUIRE ( (x==six).all() );

  x = two;
  combine<Multiply>(two,x,1,1);
  REQUIRE ( (x==four).all() );
  combine_fill_aware<Multiply>(fv,x,1,1);
  REQUIRE ( (x==four).all() );

  x = four;
  combine<Divide>(two,x,1,1);
  REQUIRE ( (x==two).all() );
  combine_fill_aware<Divide>(fv,x,1,1);
  REQUIRE ( (x==two).all() );

  x = two;
  combine<Max>(four,x,1,1);
  REQUIRE ( (x==four).all() );
  combine_fill_aware<Max>(fv,x,1,1);
  REQUIRE ( (x==four).all() );

  x = four;
  combine<Min>(two,x,1,1);
  REQUIRE ( (x==two).all() );
  x = fv_times_ten;
  combine_fill_aware<Min>(fv,x,1,1);
  REQUIRE ( (x==fv_times_ten).all() );
}

} // anonymous namespace
