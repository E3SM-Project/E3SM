#include "p3_calc_rime_density_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for computing rime density on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
