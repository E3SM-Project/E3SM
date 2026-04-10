#include "impl/zm_zm_conv_mcsp_tend_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_conv_mcsp_tend on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
