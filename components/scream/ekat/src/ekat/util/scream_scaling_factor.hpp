#ifndef SCREAM_SCALING_FACTOR_HPP
#define SCREAM_SCALING_FACTOR_HPP

#include <iostream>
#include "ekat/util/scream_rational_constant.hpp"

#include <array>

namespace scream
{

namespace util
{

struct ScalingFactor {
  const RationalConstant base;
  const RationalConstant exp;

  ScalingFactor () = delete;
  constexpr ScalingFactor (const RationalConstant& s)
   : ScalingFactor(s,1)
  {
    // Nothing to do here
  }
  constexpr ScalingFactor (const RationalConstant& b,
                           const RationalConstant& e)
   : base(std::get<0>(check_and_adjust(b,e)))
   , exp (std::get<1>(check_and_adjust(b,e)))
  {
    // Nothing to do here
  }

  constexpr ScalingFactor (const ScalingFactor&) = default;

  static constexpr ScalingFactor one () { return ScalingFactor(1); }
  static constexpr ScalingFactor zero () { return ScalingFactor(0); }

private:

  using base_and_exp = std::array<RationalConstant,2>;
  static constexpr base_and_exp
  check_and_adjust (const RationalConstant b,
                    const RationalConstant e)
  {
    // Check that we are not doing 0^0 or taking even roots of negative numbers.
    // If all good, adjust base and exp for x^0 case (b=1,e=1), and return what was requested.
    return CONSTEXPR_ASSERT( !(b==RationalConstant::zero() && e==RationalConstant::zero()) ),
           CONSTEXPR_ASSERT( !(b.num()<0 && e.den()%2==0) ),
           (e==RationalConstant::zero()
              ? base_and_exp{{RationalConstant::one(), RationalConstant::one()}}
              : base_and_exp{{b,e}} );
  }
};

constexpr bool operator== (const ScalingFactor& lhs, const ScalingFactor& rhs) {
  // Recall that, with all terms being integers,
  //    (a/b)^(c/d)==(x/y)^(w/z)
  // is equivalent to
  //    (a/b)^(cz) == (x/w)^(wd)
  return pow(lhs.base,lhs.exp.num()*rhs.exp.den())==pow(rhs.base,rhs.exp.num()*lhs.exp.den());
}

constexpr bool operator== (const ScalingFactor& lhs, const RationalConstant& rhs) {
  return lhs==ScalingFactor(rhs);
}

constexpr bool operator== (const RationalConstant& lhs, const ScalingFactor& rhs) {
  return ScalingFactor(lhs)==rhs;
}

constexpr bool operator!= (const ScalingFactor& lhs, const ScalingFactor& rhs) {
  return !(lhs==rhs);
}

constexpr ScalingFactor operator* (const ScalingFactor& lhs, const ScalingFactor& rhs) {
  // If base or exp are the same, we can use powers properties,
  // otherwise, recall that, with all terms being integers,
  //    (a/b)^(c/d) * (x/y)^(w/z)
  // is equivalent to
  //    ((a/b)^(cz) * (x/w)^(wd)) ^ (1/dz)
  return lhs.base==rhs.base 
            ? ScalingFactor(lhs.base,lhs.exp+rhs.exp)
            : (lhs.exp==rhs.exp
                ? ScalingFactor(lhs.base*rhs.base,lhs.exp)
                : ScalingFactor( pow(lhs.base,lhs.exp.num()*rhs.exp.den()) * pow(rhs.base,rhs.exp.num()*lhs.exp.den()),
                                 RationalConstant(1,lhs.exp.den()*rhs.exp.den())));
}

constexpr ScalingFactor operator* (const RationalConstant& lhs, const ScalingFactor& rhs) {
  return rhs*ScalingFactor(lhs);
}

constexpr ScalingFactor operator* (const ScalingFactor& lhs, const RationalConstant& rhs) {
  return lhs*ScalingFactor(rhs);
}

constexpr ScalingFactor operator/ (const ScalingFactor& lhs, const ScalingFactor& rhs) {
  // If base or exp are the same, we can use powers properties,
  // otherwise, recall that, with all terms being integers,
  //    (a/b)^(c/d) / (x/y)^(w/z)
  // is equivalent to
  //    ((a/b)^(cz) / (x/w)^(wd)) ^ (1/dz)
  return lhs.base==rhs.base 
            ? ScalingFactor(lhs.base,lhs.exp-rhs.exp)
            : (lhs.exp==rhs.exp
                ? ScalingFactor(lhs.base/rhs.base,lhs.exp)
                : ScalingFactor( pow(lhs.base,lhs.exp.num()*rhs.exp.den()) / pow(rhs.base,rhs.exp.num()*lhs.exp.den()),
                                 RationalConstant(1,lhs.exp.den()*rhs.exp.den())));
}

constexpr ScalingFactor operator/ (const ScalingFactor& lhs, const RationalConstant& rhs) {
  return lhs/ScalingFactor(rhs);
}

constexpr ScalingFactor operator/ (const RationalConstant& lhs, const ScalingFactor& rhs) {
  return ScalingFactor(lhs)/rhs;
}

constexpr ScalingFactor pow (const ScalingFactor& x, const int p) {
  return ScalingFactor(x.base,x.exp*p);
}

constexpr ScalingFactor pow (const ScalingFactor& x, const RationalConstant& y) {
  return ScalingFactor(x.base,x.exp*y);
}

constexpr ScalingFactor sqrt (const ScalingFactor& x) {
  return ScalingFactor(x.base,x.exp/2);
}

inline std::string to_string(const ScalingFactor& x, const Format exp_fmt = Format::Rat) {
  std::string s = to_string(x.base);
  if (x.exp!=RationalConstant::one()) {
    s +=  "^" + to_string(x.exp,exp_fmt);
  }
  return s;
}

inline std::ostream& operator<< (std::ostream& out, const ScalingFactor& s) {
  out << to_string(s);
  return out;
}


namespace prefixes {
constexpr ScalingFactor nano  = ScalingFactor(10,-9);
constexpr ScalingFactor micro = ScalingFactor(10,-6);
constexpr ScalingFactor milli = ScalingFactor(10,-3);
constexpr ScalingFactor centi = ScalingFactor(10,-2);
constexpr ScalingFactor deci  = ScalingFactor(10,-1);

constexpr ScalingFactor deca  = ScalingFactor(10, 1);
constexpr ScalingFactor hecto = ScalingFactor(10, 2);
constexpr ScalingFactor kilo  = ScalingFactor(10, 3);
constexpr ScalingFactor mega  = ScalingFactor(10, 6);
constexpr ScalingFactor giga  = ScalingFactor(10, 9);
} // namespace prefixes

} // namespace util

} // namespace scream

#endif // SCREAM_SCALING_FACTOR_HPP
