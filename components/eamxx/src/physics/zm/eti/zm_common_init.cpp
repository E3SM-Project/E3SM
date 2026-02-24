#include "impl/zm_common_init_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_common_init on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
