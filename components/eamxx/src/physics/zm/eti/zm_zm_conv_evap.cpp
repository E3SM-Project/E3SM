#include "impl/zm_zm_conv_evap_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_conv_evap on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
