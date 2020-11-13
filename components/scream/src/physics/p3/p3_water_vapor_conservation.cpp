#include "p3_water_vapor_conservation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing water_vapor_conservation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
