#ifndef SCREAM_UNITS_HPP
#define SCREAM_UNITS_HPP

#include "share/util/rational_constant.hpp"
#include "share/util/scaling_factor.hpp"

#include <sstream>

namespace scream
{

namespace units
{

// In the units namespace, it's ok to use some stuff from the util namespace
using util::RationalConstant;
using util::ScalingFactor;
using util::Format;
using namespace util::prefixes;

constexpr int NUM_BASIC_UNITS = 7;
constexpr const char* BASIC_UNITS_SYMBOLS[7] = {"m", "s", "kg", "K", "A", "mol", "cd"};

/*
 *  Units: a class to store physical units in terms of fundamental ones
 *  
 *  Units is morally storing 8 numbers:
 *   - a scaling factor
 *   - the exponents of the 7 base units
 *  So if a quantity has units kPa, it will store
 *  
 *    - a scaling factor of 1000
 *    - the exponents [ -1, -2, 1, 0, 0, 0, 0 ]
 *  
 *  since kPa = 1000 kg m^-1 s ^-2.
 */

class Units {
public:

  // No default
  constexpr Units () = delete;

  // Construct a non-dimensional quantity
  constexpr Units (const ScalingFactor& scaling)
   : m_scaling {scaling}
   , m_units{0,0,0,0,0,0,0}
   , m_exp_format (Format::Rat)
  {
    // Nothing to do here
  }

  // Construct a general quantity
  constexpr Units (const RationalConstant& lengthExp,
                   const RationalConstant& timeExp,
                   const RationalConstant& massExp,
                   const RationalConstant& temperatureExp,
                   const RationalConstant& currentExp,
                   const RationalConstant& amountExp,
                   const RationalConstant& luminousIntensityExp,
                   const ScalingFactor& scalingFactor = RationalConstant::one())
   : m_scaling {scalingFactor.base,scalingFactor.exp}
   , m_units {lengthExp,timeExp,massExp,temperatureExp,currentExp,amountExp,luminousIntensityExp}
   , m_exp_format (Format::Rat)
  {
    // Nothing to do here
  }

  constexpr Units (const Units&) = default;
  Units& operator= (const Units&) = default;

  const Units& set_exp_format (const Format fmt) {
    m_exp_format = fmt;
    return *this;
  }

private:

  friend constexpr bool operator==(const Units&,const Units&);

  friend constexpr Units operator*(const Units&,const Units&);
  friend constexpr Units operator*(const ScalingFactor&,const Units&);
  friend constexpr Units operator/(const Units&,const Units&);
  friend constexpr Units operator/(const Units&,const ScalingFactor&);
  friend constexpr Units operator/(const ScalingFactor&,const Units&);
  friend constexpr Units pow(const Units&,const RationalConstant&);
  friend constexpr Units sqrt(const Units&);
  friend constexpr Units root(const Units&,const int);

  friend std::string to_string(const Units&);

  const ScalingFactor     m_scaling;
  const RationalConstant  m_units[7];

  Format                  m_exp_format;
};

// === Operators/functions overload === //
//
// Recall: the first RationalConstant is a factor, while the remaining seven
//         are only exponents. So in u1*u2, multiply the factor, and
//         add the exponents

// --- Comparison --- //

constexpr bool operator==(const Units& lhs, const Units& rhs) {
  return lhs.m_scaling==rhs.m_scaling &&
         lhs.m_units[0]==rhs.m_units[0] &&
         lhs.m_units[1]==rhs.m_units[1] &&
         lhs.m_units[2]==rhs.m_units[2] &&
         lhs.m_units[3]==rhs.m_units[3] &&
         lhs.m_units[4]==rhs.m_units[4] &&
         lhs.m_units[5]==rhs.m_units[5] &&
         lhs.m_units[6]==rhs.m_units[6];
}

constexpr bool operator!=(const Units& lhs, const Units& rhs) {
  return !(lhs==rhs);
}

// --- Multiplicaiton --- //
constexpr Units operator*(const Units& lhs, const Units& rhs) {
  return Units(lhs.m_units[0]+rhs.m_units[0],
               lhs.m_units[1]+rhs.m_units[1],
               lhs.m_units[2]+rhs.m_units[2],
               lhs.m_units[3]+rhs.m_units[3],
               lhs.m_units[4]+rhs.m_units[4],
               lhs.m_units[5]+rhs.m_units[5],
               lhs.m_units[6]+rhs.m_units[6],
               lhs.m_scaling*rhs.m_scaling);
}
constexpr Units operator*(const ScalingFactor& lhs, const Units& rhs) {
  return Units(rhs.m_units[0],
               rhs.m_units[1],
               rhs.m_units[2],
               rhs.m_units[3],
               rhs.m_units[4],
               rhs.m_units[5],
               rhs.m_units[6],
               lhs*rhs.m_scaling);
}
constexpr Units operator*(const Units& lhs, const ScalingFactor& rhs) {
  return rhs*lhs;
}
constexpr Units operator*(const RationalConstant& lhs, const Units& rhs) {
  return ScalingFactor(lhs)*rhs;
}
constexpr Units operator*(const Units& lhs, const RationalConstant& rhs) {
  return lhs*ScalingFactor(rhs);
}

// --- Division --- //
constexpr Units operator/(const Units& lhs, const Units& rhs) {
  return Units(lhs.m_units[0]-rhs.m_units[0],
               lhs.m_units[1]-rhs.m_units[1],
               lhs.m_units[2]-rhs.m_units[2],
               lhs.m_units[3]-rhs.m_units[3],
               lhs.m_units[4]-rhs.m_units[4],
               lhs.m_units[5]-rhs.m_units[5],
               lhs.m_units[6]-rhs.m_units[6],
               lhs.m_scaling/rhs.m_scaling);
}
constexpr Units operator/(const Units& lhs, const ScalingFactor& rhs) {
  return Units(lhs.m_units[0],
               lhs.m_units[1],
               lhs.m_units[2],
               lhs.m_units[3],
               lhs.m_units[4],
               lhs.m_units[5],
               lhs.m_units[6],
               lhs.m_scaling/rhs);
}
constexpr Units operator/(const ScalingFactor& lhs, const Units& rhs) {
  return Units(1-rhs.m_units[0],
               1-rhs.m_units[1],
               1-rhs.m_units[2],
               1-rhs.m_units[3],
               1-rhs.m_units[4],
               1-rhs.m_units[5],
               1-rhs.m_units[6],
               lhs/rhs.m_scaling);
}
constexpr Units operator/(const RationalConstant& lhs, const Units& rhs) {
  return ScalingFactor(lhs)/rhs;
}
constexpr Units operator/(const Units& lhs, const RationalConstant& rhs) {
  return lhs/ScalingFactor(rhs);
}

// --- Powers and roots --- //

constexpr Units pow(const Units& x, const RationalConstant& p) {
  return Units(x.m_units[0]*p,
               x.m_units[1]*p,
               x.m_units[2]*p,
               x.m_units[3]*p,
               x.m_units[4]*p,
               x.m_units[5]*p,
               x.m_units[6]*p,
               pow(x.m_scaling,p));
}

constexpr Units sqrt(const Units& x) {
  return Units(x.m_units[0] / 2,
               x.m_units[1] / 2,
               x.m_units[2] / 2,
               x.m_units[3] / 2,
               x.m_units[4] / 2,
               x.m_units[5] / 2,
               x.m_units[6] / 2,
               pow(x.m_scaling,RationalConstant(1,2)));
}

inline std::string to_string(const Units& x) {
  std::string s;
  int num_non_trivial = 0;
  for (int i=0; i<NUM_BASIC_UNITS; ++i) {
    if (x.m_units[i].num()==0) {
      continue;
    }
    ++num_non_trivial;
    s += " ";
    s += BASIC_UNITS_SYMBOLS[i];
    if (x.m_units[i]!=RationalConstant::one()) {
      s += "^" + to_string(x.m_units[i],x.m_exp_format);
    }
  }

  // Prepend the scaling only if not one or dimensionless unit
  if (x.m_scaling!=RationalConstant::one() || num_non_trivial==0) {
    s = to_string(x.m_scaling,x.m_exp_format) + s;
  }
  return s;
}

inline std::ostream& operator<< (std::ostream& out, const Units& x) {
  out << to_string(x);
  return out;
}

// ================== SHORT NAMES FOR COMMON PHYSICAL UNITS =============== //

// Note to developers:
// I added the 'most common' units. I avoided 'Siemes', since
// the symbol is 'S', which is way too close to 's' (seconds).
// TODO: should we add 'common' scaled ones, such as km=kilo*m?

// === FUNDAMENTAL === //

constexpr Units m   = Units(1,0,0,0,0,0,0);
constexpr Units s   = Units(0,1,0,0,0,0,0);
constexpr Units kg  = Units(0,0,1,0,0,0,0);
constexpr Units K   = Units(0,0,0,1,0,0,0);
constexpr Units A   = Units(0,0,0,0,1,0,0);
constexpr Units mol = Units(0,0,0,0,0,1,0);
constexpr Units cd  = Units(0,0,0,0,0,0,1);

// === DERIVED === //

// Thermomechanics
constexpr auto day  = 86400*s;      // day          (time)
constexpr auto year = 365*day;      // year         (time)
constexpr auto g    = milli*kg;     // gram         (mass)
constexpr auto N    = kg*m/(s*s);   // newton       (force)
constexpr auto dyn  = N/(10000);    // dyne         (force)
constexpr auto Pa   = N/(m*m);      // pascal       (pressure)
constexpr auto bar  = 10000*Pa;     // bar          (pressure)
constexpr auto atm  = 101325*Pa;    // atmosphere   (pressure)
constexpr auto J    = N*m;          // joule        (energy)
constexpr auto W    = J/s;          // watt         (power)

// Electro-magnetism
constexpr auto C    = A*s;          // coulomb      (charge)
constexpr auto V    = J/C;          // volt         (voltage)
constexpr auto T    = N/(A*m);      // tesla        (magnetic field)
constexpr auto F    = C/V;          // farad        (capacitance)
constexpr auto Wb   = V*s;          // weber        (magnetic flux)
constexpr auto H    = Wb/A;         // henri        (inductance)
constexpr auto Sv   = J/kg;         // sievert      (radiation dose)
constexpr auto rem  = Sv/100;       // rem          (radiation dose)
constexpr auto Hz   = 1/s;          // hertz        (frequency)

} // namespace units

} // namespace scream

#endif // SCREAM_UNITS_HPP
