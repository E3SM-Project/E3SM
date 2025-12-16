#ifndef EAMXX_QUANTITY_HPP
#define EAMXX_QUANTITY_HPP

#include <ekat_math_utils.hpp>
#include <ekat_units.hpp>

namespace scream {

template<typename Scalar>
struct Quantity
{
  using Units = ekat::units::Units;

  constexpr Quantity() : Quantity(ekat::invalid<Scalar>()) {}
  constexpr Quantity(const Scalar v) : value(v), units(Units::invalid()) {}
  constexpr Quantity(const Scalar v, const Units& u) : value(v), units(u) {}

  Scalar  value;
  Units   units;
};

// Nano-util, that allows one to not having to bother to figure out the exact type
// of the scalar when creating a quantity object
template<typename S>
constexpr Quantity<S> create_quantity(const S s, const ekat::units::Units& u)
{
  return Quantity<S>(s,u);
}
// ---------------- Operators overloads ------------------ //

template<typename S1, typename S2>
constexpr auto operator* (const Quantity<S1>& l, const Quantity<S2>& r)
{
  return create_quantity(l.value*r.value,l.units*r.units);
}
template<typename S1, typename S2>
constexpr auto operator* (const S1& l, const Quantity<S1>& r)
{
  return create_quantity(l*r.value,r.units);
}
template<typename S1, typename S2>
constexpr auto operator* (const Quantity<S1>& l, const S2 r)
{
  return create_quantity(l.value*r,l.units);
}
template<typename S1, typename S2>
constexpr auto operator/ (const Quantity<S1>& l, const Quantity<S2>& r)
{
  return create_quantity(l.value/r.value,l.units*r.units);
}
template<typename S1, typename S2>
constexpr auto operator/ (const S1 l, const Quantity<S2>& r)
{
  return create_quantity(l/r.value,r.units);
}
template<typename S1, typename S2>
constexpr auto operator/ (const Quantity<S1>& l, const S2 r) {
  return create_quantity(l.value/r,l.units);
}

} // namespace scream

#endif // EAMXX_QUANTITY_HPP
