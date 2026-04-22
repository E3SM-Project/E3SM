#include "impl/zm_cloud_properties_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_cloud_properties on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
