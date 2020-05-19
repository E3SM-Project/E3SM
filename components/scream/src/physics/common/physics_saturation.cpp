#include "p3_functions_math_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace physics {

/*
 * Explicit instantiation for doing p3 math functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace physics
} // namespace scream
