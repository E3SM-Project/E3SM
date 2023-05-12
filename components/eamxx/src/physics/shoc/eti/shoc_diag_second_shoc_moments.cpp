#include "shoc_diag_second_shoc_moments_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing diag_second_shoc_moments on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
