#ifndef SCREAM_RATIONAL_CONSTANT_HPP
#define SCREAM_RATIONAL_CONSTANT_HPP

#include <stdexcept>
#include <sstream>

#include "ekat/scream_assert.hpp"
#include "ekat/scream_types.hpp"

namespace scream
{

namespace util
{

enum class Format {
  Float,
  Rat,
  Auto
};

inline constexpr Format operator&& (const Format fmt1, const Format fmt2) {
  return fmt1==Format::Rat || fmt2==Format::Rat ? Format::Rat
                    : (fmt1==Format::Auto ? fmt2 : fmt1);
}

constexpr long long abs (const long long a) {
  return a>=0 ? a : -a;
}

constexpr long long gcd (const long long a, const long long b) {
  return b==0
          ? a
          : gcd(b, a % b);
}

// Elementary representation of a rational number
// Note: we do not allow construction from a double (even if
//       the double is a simple fraction, like 1.5), since
//       its implementation (with correct checks and safeguards)
//       would probably require too much effort, especially
//       considered the scope of the usage of this class

class RationalConstant {
public:

  // No default
  RationalConstant () = delete;

  // No construction from Real
  RationalConstant (const Real x) = delete;

  // Construction from long, means long/1
  template<typename IntType>
  constexpr RationalConstant (const IntType n,
                              typename std::enable_if<std::is_integral<IntType>::value,Format>::type fmt = Format::Float)
   : RationalConstant (n,1,fmt)
  {
    // Nothing to do here
  }

  template<typename IntType1, typename IntType2>
  constexpr RationalConstant (const IntType1 n, const IntType2 d,
                              typename std::enable_if<std::is_integral<IntType1>::value &&
                                                      std::is_integral<IntType1>::value,Format>::type fmt = Format::Float)
   : m_num(fix_num(n,d))
   , m_den(fix_den(n,d))
   , m_fmt(fmt)
  {
    // Nothing to do here
  }
  constexpr RationalConstant (const RationalConstant&) = default;

  constexpr RationalConstant operator- () const {
    return RationalConstant(-m_num,m_den,m_fmt);
  }

  constexpr long long num () const { return m_num; }
  constexpr long long den () const { return m_den; }

  const RationalConstant& set_format (const Format fmt) {
    m_fmt = fmt;
    return *this;
  }

  constexpr Format get_format () const { return m_fmt; }

  static constexpr RationalConstant one () { return RationalConstant(1); }
  static constexpr RationalConstant zero () { return RationalConstant(0); }

private:

  // These two are used to help reduce a/b to lowest terms
  static constexpr long long fix_num(const long long n, const long long d) {
    return CONSTEXPR_ASSERT(d!=0),
             n==0 ? n : n / gcd(abs(n),abs(d));
  }

  static constexpr long long fix_den(const long long n, const long long d) {
    return CONSTEXPR_ASSERT(d!=0),
            d<0 ? fix_den(n,-d) : d / gcd(abs(n),abs(d));
  }

  const long long m_num;
  const long long m_den;
  Format m_fmt;
};

constexpr bool operator== (const RationalConstant& lhs, const RationalConstant& rhs) {
  // Cross multiply, in case somehow (should not be possible, unless
  // someone changed the implementation) someone managed to get lhs and/or
  // rhs to not be already reduced to lower terms
  return lhs.num()*rhs.den()==lhs.den()*rhs.num();
}

constexpr bool operator!= (const RationalConstant& lhs, const RationalConstant& rhs) {
  return !(lhs==rhs);
}

constexpr RationalConstant operator+ (const RationalConstant& lhs, const RationalConstant& rhs) {
  return RationalConstant (lhs.num()*rhs.den() + lhs.den()*rhs.num(),
                           lhs.den()*rhs.den(),
                           lhs.get_format()&&rhs.get_format());
}

constexpr RationalConstant operator- (const RationalConstant& lhs, const RationalConstant& rhs) {
  return RationalConstant (lhs.num()*rhs.den() - lhs.den()*rhs.num(),
                           lhs.den()*rhs.den(),
                           lhs.get_format()&&rhs.get_format());
}

constexpr RationalConstant operator* (const RationalConstant& lhs, const RationalConstant& rhs) {
  return RationalConstant (lhs.num()*rhs.num(),
                           lhs.den()*rhs.den(),
                           lhs.get_format()&&rhs.get_format());
}

constexpr RationalConstant operator/ (const RationalConstant& lhs, const RationalConstant& rhs) {
  return RationalConstant (lhs.num()*rhs.den(),
                           lhs.den()*rhs.num(),
                           lhs.get_format()&&rhs.get_format());
}

// WARNING! If exp>>1, this generates a huge recursion. Use carefully!!
template<typename IntType>
inline constexpr RationalConstant pow (const typename std::enable_if<std::is_integral<IntType>::value,RationalConstant>::type& x, const IntType p) {
  // Three main cases:
  //  - p<0: compute the -p power of 1/x
  //  - p=0: base case, return 1 if x!=0, throw if x==0
  //  - recursion step: x^p = x * x^{p-1}
  // Note: recall that expressions like "blah1, blah2" executes both blah1 and blah2 and return blah2.
  return p<0 ? pow(1/x,-p)
             : (p==0 ? CONSTEXPR_ASSERT(x.num()!=0), RationalConstant::one()
                     : ( (p&1)!=0 ? x*pow(x*x,p>>1) : pow(x*x,p>>1)));
}

inline std::string to_string (const RationalConstant& rat, const Format fmt = Format::Auto) {
  // Note: using std::to_string(double) causes insertion of trailing zeros,
  //       and the insertion of decimal part for integers.
  //       E.g., to_string(1.0/2) leads to 0.5000000
  std::stringstream ss;
  switch (fmt) {
    case Format::Float:
      ss << static_cast<Real>(rat.num())/rat.den();
      break;
    case Format::Rat:
      ss << rat.num();
      if (rat.den()!=1) {
        ss << "/" << rat.den();
      }
      break;
    case Format::Auto:
      scream_require_msg(rat.get_format()!=Format::Auto, "Error! No format specified in the rational constant.\n");
      ss << to_string(rat,rat.get_format());
      break;
    default:
      scream_require_msg(false,"Error! Unrecognized format for printing RationalConstant.\n");
      
  }
  return ss.str();
}

inline std::ostream& operator<< (std::ostream& out, const RationalConstant& rat) {
  out << to_string(rat);
  return out;
}

} // namespace util

} // namespace scream

#endif // SCREAM_RATIONAL_CONSTANT_HPP
