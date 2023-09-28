#include "p3_ice_sed_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing ice sedimentation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
