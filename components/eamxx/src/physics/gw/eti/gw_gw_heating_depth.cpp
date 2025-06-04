#include "impl/gw_gw_heating_depth_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_heating_depth on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
