#include "shoc_isotropic_ts_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing isotropic_ts on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
