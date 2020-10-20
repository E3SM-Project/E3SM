#include "shoc_update_prognostics_implicit_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing update_prognostics_implicit on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
