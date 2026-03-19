#include "impl/gw_convect_init_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_convect_init on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
