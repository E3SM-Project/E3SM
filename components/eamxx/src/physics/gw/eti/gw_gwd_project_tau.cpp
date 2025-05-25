#include "impl/gw_gwd_project_tau_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing gwd_project_tau on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
