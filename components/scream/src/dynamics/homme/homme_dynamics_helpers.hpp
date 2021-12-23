#ifndef SCREAM_HOMME_DYNAMICS_HELPERS_HPP
#define SCREAM_HOMME_DYNAMICS_HELPERS_HPP

#include "ekat/ekat_scalar_traits.hpp"

// Homme includes
#include "Types.hpp"

namespace ekat
{

// Specialize ScalarProperties on Homme's Scalar type
template<>
struct ScalarTraits<Homme::Scalar> {
  using value_type  = Homme::Scalar;
  using scalar_type = Homme::Real;
  using inner_traits = ScalarTraits<Homme::Real>;

  static constexpr bool is_simd = true;

  static std::string name () {
    return "Homme::Scalar";
  }

  KOKKOS_INLINE_FUNCTION
  static const value_type quiet_NaN () {
    return value_type(inner_traits::quiet_NaN());
  }

  KOKKOS_INLINE_FUNCTION
  static const value_type invalid () {
    return value_type(inner_traits::invalid());
  }
};

} // namespace ekat

#endif // SCREAM_HOMME_DYNAMICS_HELPERS_HPP
