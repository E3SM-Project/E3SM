#include "compute_tms_impl.hpp"

namespace scream {
namespace tms {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace tms
} // namespace scream
