#include "impl/zm_conv_main_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_conv_main on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
