#include "shoc_linear_interp_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing linear_interp on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
