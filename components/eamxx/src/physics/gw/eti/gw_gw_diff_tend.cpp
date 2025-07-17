#include "impl/gw_gw_diff_tend_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_diff_tend on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
