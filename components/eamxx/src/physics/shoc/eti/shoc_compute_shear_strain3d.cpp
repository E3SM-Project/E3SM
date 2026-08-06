#include "shoc_compute_shear_strain3d_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing compute_shear_strain3d on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
