#include "p3_table3_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing table functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
