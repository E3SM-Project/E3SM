#include "p3_update_prognostics_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing update prognostics functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
