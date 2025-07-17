#include "impl/gw_gw_beres_src_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_beres_src on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
