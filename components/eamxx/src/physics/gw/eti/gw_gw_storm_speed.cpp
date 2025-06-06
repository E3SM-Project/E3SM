#include "impl/gw_gw_storm_speed_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_storm_speed on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
