#include "impl/gw_gwd_precalc_rhoi_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gwd_precalc_rhoi on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
