#include "impl/gw_gwd_compute_tendencies_from_stress_divergence_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gwd_compute_tendencies_from_stress_divergence on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
