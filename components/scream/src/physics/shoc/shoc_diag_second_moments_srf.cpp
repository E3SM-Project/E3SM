#include "shoc_diag_second_moments_srf_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

