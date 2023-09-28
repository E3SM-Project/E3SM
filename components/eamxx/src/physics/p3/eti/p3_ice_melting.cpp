#include "p3_ice_melting_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing ice melting on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
