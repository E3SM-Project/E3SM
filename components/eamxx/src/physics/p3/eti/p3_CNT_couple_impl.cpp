#include "p3_CNT_couple_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing conservation functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
