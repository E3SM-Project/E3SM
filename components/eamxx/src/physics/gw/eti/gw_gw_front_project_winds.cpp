#include "impl/gw_gw_front_project_winds_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_front_project_winds on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
