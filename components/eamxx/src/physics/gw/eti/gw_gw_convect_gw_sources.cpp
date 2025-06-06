#include "impl/gw_gw_convect_gw_sources_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_convect_gw_sources on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
