#include "p3_ice_supersat_conservation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing ice_supersat_conservation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
