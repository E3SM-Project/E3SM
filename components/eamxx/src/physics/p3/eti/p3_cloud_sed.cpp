#include "p3_cloud_sed_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instatiation for doing cloud sedimentation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
