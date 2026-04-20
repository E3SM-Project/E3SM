#include "impl/zm_conv_mcsp_calculate_shear_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_conv_mcsp_calculate_shear on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
