#include "p3_main_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
