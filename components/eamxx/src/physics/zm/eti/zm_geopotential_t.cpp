#include "impl/zm_geopotential_t_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_geopotential_t on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
