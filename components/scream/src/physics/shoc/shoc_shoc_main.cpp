#include "shoc_shoc_main_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing shoc_main on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
