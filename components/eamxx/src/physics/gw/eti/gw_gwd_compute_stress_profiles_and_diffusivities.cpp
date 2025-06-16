#include "impl/gw_gwd_compute_stress_profiles_and_diffusivities_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gwd_compute_stress_profiles_and_diffusivities on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
