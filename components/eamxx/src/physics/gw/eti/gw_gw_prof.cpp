#include "impl/gw_gw_prof_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_prof on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
