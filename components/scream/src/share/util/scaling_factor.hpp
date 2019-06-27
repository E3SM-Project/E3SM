#ifndef SCREAM_SCALING_FACTOR_HPP
#define SCREAM_SCALING_FACTOR_HPP

#include <iostream>
#include "share/util/rational_constant.hpp"

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
   : base(b)
   , exp (e)
  {
    // Nothing to do here
  }
  constexpr ScalingFactor (const ScalingFactor&) = default;
  ScalingFactor& operator= (const ScalingFactor&) = default;
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
// auto l = pow(lhs.base,lhs.exp.num()*rhs.exp.den());
// std::cout << "l = " << l.set_format(Rat) << "\n";
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
constexpr ScalingFactor kilo  = ScalingFactor(10, 3);
constexpr ScalingFactor mega  = ScalingFactor(10, 6);
constexpr ScalingFactor giga  = ScalingFactor(10, 9);
} // namespace prefixes

} // namespace util

} // namespace scream

#endif // SCREAM_SCALING_FACTOR_HPP
