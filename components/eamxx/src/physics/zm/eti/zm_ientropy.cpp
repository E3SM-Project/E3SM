#include "impl/zm_ientropy_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing ientropy on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
