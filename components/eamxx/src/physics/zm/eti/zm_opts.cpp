#include "impl/zm_opts_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_opts_init on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
