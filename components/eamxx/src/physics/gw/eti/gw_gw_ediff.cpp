#include "impl/gw_gw_ediff_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gw_ediff on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
