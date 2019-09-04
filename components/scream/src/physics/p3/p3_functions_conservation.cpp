#include "p3_functions_conservation_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instatiation for doing p3 conservation functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream