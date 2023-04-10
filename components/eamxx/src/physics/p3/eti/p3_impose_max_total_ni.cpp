#include "p3_impose_max_total_ni_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing impose_max_total_ni_impl on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
