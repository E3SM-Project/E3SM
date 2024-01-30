#include "impl/dp_crm_resolved_turb_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing crm_resolved_turb on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
