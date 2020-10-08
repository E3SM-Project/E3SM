#include "shoc_diag_second_moments_lbycond_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing diag_second_moments_lbycond on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
