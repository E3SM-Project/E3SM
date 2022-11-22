#include "p3_dsd2_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instatiation for doing dsd2 functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
