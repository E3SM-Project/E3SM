#include "impl/zm_entropy_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing entropy on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
