#ifndef SCREAM_HOMME_DYNAMICS_HELPERS_HPP
#define SCREAM_HOMME_DYNAMICS_HELPERS_HPP

// Scream includes
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"

// Homme includes
#include "Types.hpp"

namespace scream
{

namespace util
{

// Specialize ScalarProperties on Homme's Scalar type
template<>
struct ScalarProperties<Homme::Scalar> {
  using scalar_type = Homme::Real;
  static constexpr bool is_pack = true;
};

template<>
struct TypeName<Homme::Scalar> {
  static std::string name () {
    return "Homme::Scalar";
  }
};


} // namespace util

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HELPERS_HPP
