#include "impl/gw_gw_common_init_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_common_init on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
