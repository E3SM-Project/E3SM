#include "impl/zm_find_mse_max_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing find_mse_max on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
