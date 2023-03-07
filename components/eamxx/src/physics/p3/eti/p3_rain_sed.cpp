#include "p3_rain_sed_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing rain sedimentation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
