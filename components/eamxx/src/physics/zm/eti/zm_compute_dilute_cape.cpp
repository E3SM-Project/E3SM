#include "impl/zm_compute_dilute_cape_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing compute_dilute_cape on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
