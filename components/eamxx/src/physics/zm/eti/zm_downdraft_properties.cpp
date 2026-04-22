#include "impl/zm_downdraft_properties_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_downdraft_properties on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
