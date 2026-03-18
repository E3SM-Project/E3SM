#include "impl/zm_invert_entropy_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing invert_entropy on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
