#include "shoc_diag_second_moments_ubycond_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 *  * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

