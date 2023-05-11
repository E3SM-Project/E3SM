#include "p3_autoconversion_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing autoconversion on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
