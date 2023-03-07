#include "shoc_calc_shoc_vertflux_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing vertflux on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
