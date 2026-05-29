#include "impl/zm_state_update_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_state_update on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
