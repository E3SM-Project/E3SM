#include "impl/zm_closure_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_closure on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
