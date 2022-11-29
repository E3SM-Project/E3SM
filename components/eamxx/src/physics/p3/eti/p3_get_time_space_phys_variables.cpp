#include "p3_get_time_space_phys_variables_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing get_time_space_phys_variables function on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
