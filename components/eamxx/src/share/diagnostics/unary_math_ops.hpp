#ifndef EAMXX_UNARY_MATH_OPS_HPP
#define EAMXX_UNARY_MATH_OPS_HPP

#include <Kokkos_Core.hpp>

namespace scream {

namespace math_ops {

struct Log {
  static constexpr const char* name() { return "log"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::log(x); }
};

struct Exp {
  static constexpr const char* name() { return "exp"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::exp(x); }
};

struct Sin {
  static constexpr const char* name() { return "sin"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::sin(x); }
};

struct Cos {
  static constexpr const char* name() { return "cos"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::cos(x); }
};

struct Tan {
  static constexpr const char* name() { return "tan"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::tan(x); }
};

struct Atan {
  static constexpr const char* name() { return "atan"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::atan(x); }
};

struct Sqrt {
  static constexpr const char* name() { return "sqrt"; }
  KOKKOS_INLINE_FUNCTION double operator()(const double x) const { return Kokkos::sqrt(x); }
};

} // namespace math_ops

} // namespace scream

#endif // EAMXX_UNARY_MATH_OPS_HPP
