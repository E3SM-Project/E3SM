#include "rayleigh_friction_impl.hpp"

namespace scream {
namespace rayleigh_friction {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace rayleigh_friction
} // namespace scream
