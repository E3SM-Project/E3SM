#include "shoc_compute_shoc_vapor_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing compute_shoc_vapor on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
