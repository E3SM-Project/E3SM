#include "p3_nc_conservation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing nc_conservation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
