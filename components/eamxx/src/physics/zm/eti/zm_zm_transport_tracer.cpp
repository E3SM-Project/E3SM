#include "impl/zm_zm_transport_tracer_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_transport_tracer on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
