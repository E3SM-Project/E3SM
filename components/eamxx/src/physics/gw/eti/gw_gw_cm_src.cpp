#include "impl/gw_gw_cm_src_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_cm_src on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
