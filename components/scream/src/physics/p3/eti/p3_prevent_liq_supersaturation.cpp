#include "p3_prevent_liq_supersaturation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing prevent_liq_supersaturation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
