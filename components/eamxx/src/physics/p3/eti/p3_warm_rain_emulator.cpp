#include "p3_warm_rain_emulator_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for merging emulated warm-rain rates on Reals using
 * the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
