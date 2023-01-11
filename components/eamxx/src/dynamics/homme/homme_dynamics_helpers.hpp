#ifndef SCREAM_HOMME_DYNAMICS_HELPERS_HPP
#define SCREAM_HOMME_DYNAMICS_HELPERS_HPP

#include "Context.hpp"

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

namespace scream {

// Since we have multiple access points to Homme's context,
// Use this tiny helper class to allow correct finalization
// of the context. Last one to exit shuts the door.
struct HommeContextUser {
  void add_user () { ++counter; }
  void remove_user () {
    --counter;
    if (counter==0) {
      Homme::Context::singleton().finalize_singleton();
    }
  }
  int get_counter () const { return counter; }
  static HommeContextUser& singleton() {
    static HommeContextUser hcu;
    return hcu;
  }
private:
  HommeContextUser () : counter(0) {}
  int counter;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HELPERS_HPP
