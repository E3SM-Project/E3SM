#ifndef SCREAM_HOMME_DYNAMICS_HELPERS_HPP
#define SCREAM_HOMME_DYNAMICS_HELPERS_HPP

#include "Context.hpp"

#include <ekat_scalar_traits.hpp>

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

  static HommeContextUser& singleton() {
    static HommeContextUser hcu;
    return hcu;
  }
private:
  HommeContextUser () = default;
  int counter = 0;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HELPERS_HPP
